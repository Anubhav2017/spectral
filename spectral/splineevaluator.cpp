#include "splineevaluator.h"

#include "parameteroptimizer.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <Dense>

SplineEvaluator::SplineEvaluator(BSplinePatchNetwork *bSplinePatchNetwork, Eigen::SparseMatrix<double> *energyMatrix, int numberOfEvalPointsForG1, bool doOptimizationForDataError):
    m_bSplinePatchNetwork(bSplinePatchNetwork), m_cellMesh(bSplinePatchNetwork->getCellMesh()), m_mesh(m_cellMesh->getOriginalMesh()), m_energyMatrix(energyMatrix)
{
    double minX, maxX, minY, maxY, minZ, maxZ;
    m_mesh->calculateMinMaxVertexPositions(minX, maxX, minY, maxY, minZ, maxZ);

    m_meshSizeX = maxX - minX;
    m_meshSizeY = maxY - minY;
    m_meshSizeZ = maxZ - minZ;

    m_mesh->calculateMinimumBoundingSphere(m_boundingSphereRadius, m_boundingSphereCenter);

    this->evaluateG1Errors(numberOfEvalPointsForG1);
    this->evaluateLSErrors(doOptimizationForDataError);
    this->evaluatePatchEnergyValues();
    this->evaluatePatchAreas();
}

void SplineEvaluator::evaluateLSErrors(bool doOptimizationForDataError)
{
    const int numberOfCellFaces = m_cellMesh->getNumberOfFaces();
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    const int numberOfMeshVertices = m_mesh->getNumberOfVertices();

    m_dataErrors = QVector<double>(numberOfMeshVertices, 999999);
    m_projectedCoordinates = QVector<QVector3D>(numberOfMeshVertices);

    //calculate LS error
    int totalNumberOfDatapoints = 0;
    m_totalSquaredError = 0;
    for (int iCell = 0; iCell < numberOfCellFaces; iCell++) {
        CellFace *cf = m_cellMesh->getFace(iCell);
        Parameterization *param = m_cellMesh->getQuadCellParameterization(iCell);
        const int numberOfDatapoints = cf->getMeshVertices()->size();
        totalNumberOfDatapoints += numberOfDatapoints;

        #pragma omp parallel for
        for (int k = 0; k < numberOfDatapoints; k++) {
            const QVector3D vertexPos = cf->getMeshVertices()->at(k)->getPosition();
            const QVector2D paramValue = param->getParameter(cf->getMeshVertices()->at(k));
            QVector3D bsplinePos = m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluate(paramValue.x(), paramValue.y());
            QVector3D errorVector =  bsplinePos - vertexPos;
            double dist = errorVector.length();
            QVector3D normal = QVector3D::crossProduct(m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluateTu(paramValue.x(), paramValue.y()),
                                                       m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluateTv(paramValue.x(), paramValue.y()));

            if (doOptimizationForDataError) {
                //do orthogonal projection to find a better parameter value
                const QVector2D newParamValue = ParameterOptimizer::optimizeParameter_surfacePoint_projection_unbounded(m_bSplinePatchNetwork->getBSplineSurface(iCell), vertexPos, paramValue);

                //check if parameter value is still in [0,1]x[0,1]
                const bool overEdge[4] = {newParamValue.y() < 0, newParamValue.x() > 1, newParamValue.y() > 1, newParamValue.x() < 0};

                int internalCvId = -1;
                int internalCeId = -1;
                for (int iEdge = 0; iEdge < 4; iEdge++) {
                    const int iEdgeNext = (iEdge + 1) % 4;
                    if (overEdge[iEdge]) {
                        internalCeId = iEdge;
                        if (overEdge[iEdgeNext])
                            internalCvId = iEdgeNext;   //cv (i+1) is between ce i and ce (i+1)
                    }
                }

                //if optimized parameter is outside, check if one of the neighboring patches might offer a better fit
                QVector<int> faceIDsToConsider;
                if (internalCvId != -1)
                    faceIDsToConsider = *(cf->getVertex(internalCvId)->getIncidentFaceIds());
                else if (internalCeId != -1)
                    faceIDsToConsider << cf->getEdge(internalCeId)->getIncidentFaceId(0) << cf->getEdge(internalCeId)->getIncidentFaceId(1);

                foreach (int incidentCellFaceId, faceIDsToConsider) {   //includes also the original face (but with nlo we might find a better result)
                    const QVector2D optParam = ParameterOptimizer::optimizeParameter_surfacePoint_nlo_bounded(m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceId), vertexPos);
                    const QVector3D optBSplinePos = m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceId)->evaluate(optParam.x(), optParam.y());
                    const QVector3D optErrorVector = optBSplinePos - vertexPos;
                    const double optDist = optErrorVector.length();
                    if (optDist < dist) {
                        bsplinePos = optBSplinePos;
                        errorVector = optErrorVector;
                        dist = optDist;
                        normal = QVector3D::crossProduct(m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceId)->evaluateTu(optParam.x(), optParam.y()),
                                                         m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceId)->evaluateTv(optParam.x(), optParam.y()));
                    }
                }

                if (!overEdge[0] && !overEdge[1] && !overEdge[2] && !overEdge[3]) {
                    const QVector3D newBSplinePos = m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluate(newParamValue.x(), newParamValue.y());
                    const QVector3D newErrorVector =  newBSplinePos - vertexPos;
                    const double newDist = newErrorVector.length();
                    if (newDist < dist) {
                        bsplinePos = newBSplinePos;
                        errorVector = newErrorVector;
                        dist = newDist;
                        normal = QVector3D::crossProduct(m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluateTu(newParamValue.x(), newParamValue.y()),
                                                         m_bSplinePatchNetwork->getBSplineSurface(iCell)->evaluateTv(newParamValue.x(), newParamValue.y()));
                    }
                }
            }

            if (QVector3D::dotProduct(normal, errorVector) < 0)
                    dist *= -1;

            const int vertexId = cf->getMeshVertices()->at(k)->getId();
            m_dataErrors[vertexId] = dist;
            m_projectedCoordinates[vertexId] = bsplinePos;
        }
    }

    //TODO find optimal parameters for cell edge plvs
    for (int iCE = 0; iCE < numberOfCellEdges; iCE++) {
        CellEdge *ce = m_cellMesh->getEdge(iCE);
        const int cfId = ce->getIncidentFaceId(0);
        Parameterization *param = m_cellMesh->getQuadCellParameterization(cfId);
        const int numPLV = ce->getPolylineVertices()->size();
        for (int i = 0; i < numPLV; i++) {
            PolylineVertex *plv = ce->getPolylineVertices()->at(i);
            if (plv->isMeshVertex()) {
                totalNumberOfDatapoints++;
                const QVector2D paramValue = param->getParameter(plv);
                const QVector3D bsplinePos = m_bSplinePatchNetwork->getBSplineSurface(cfId)->evaluate(paramValue.x(), paramValue.y());
                const QVector3D errorVector = bsplinePos - plv->getPosition();
                double dist = errorVector.length();
                QVector3D normal = QVector3D::crossProduct(m_bSplinePatchNetwork->getBSplineSurface(cfId)->evaluateTu(paramValue.x(), paramValue.y()),
                                                           m_bSplinePatchNetwork->getBSplineSurface(cfId)->evaluateTv(paramValue.x(), paramValue.y()));
                if (QVector3D::dotProduct(normal, errorVector) < 0)
                        dist *= -1;
                const int vertexId = plv->getMeshVertex()->getId();
                m_dataErrors[vertexId] = dist;
                m_projectedCoordinates[vertexId] = bsplinePos;
            }
        }
    }

    m_maxDataError = -1;
    for (int k = 0; k < numberOfMeshVertices; k++) {
        double dist = fabs(m_dataErrors[k]);
        m_totalSquaredError += dist * dist;

        if (dist > m_maxDataError)
            m_maxDataError = dist;
    }

    if (numberOfMeshVertices != totalNumberOfDatapoints)
        qDebug() << "Something went wrong with evaluating the data error. Number of mesh vertices:" << numberOfMeshVertices << "Number of evaluated points:" << totalNumberOfDatapoints;

    m_rmsError = sqrt(m_totalSquaredError/numberOfMeshVertices);

    const int numberOfMeshEdges = m_mesh->getNumberOfEdges();
    m_averageEdgeLength = 0;
    for (int i = 0; i < numberOfMeshEdges; i++)
        m_averageEdgeLength += m_mesh->getEdge(i)->getDistanceBetweenVertices();

    m_averageEdgeLength /= numberOfMeshEdges;

    double maxSize = m_meshSizeX;
    if (m_meshSizeY > maxSize)
        maxSize = m_meshSizeY;
    if (m_meshSizeZ > maxSize)
        maxSize = m_meshSizeZ;

    m_boundingBoxScaledRmsError = m_rmsError/maxSize;
    m_boundingBoxScaledMaxDataError = m_maxDataError/maxSize;

    m_boundingSphereScaledRmsErrorPercent = m_rmsError/(m_boundingSphereRadius * 2) * 100;  //scaled to have diameter 1
    m_boundingSphereScaledMaxDataErrorPercent = m_maxDataError/(m_boundingSphereRadius * 2) * 100;
}

void SplineEvaluator::evaluateG1Errors(int numberOfG1SamplePoints)
{
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();

    m_individualDetErrors = QVector<QVector<double> >(numberOfCellEdges, QVector<double>(numberOfG1SamplePoints));
    m_integratedDetErrors = QVector<double>(numberOfCellEdges, 0);
    m_individualNormalErrors = QVector<QVector<double> >(numberOfCellEdges, QVector<double>(numberOfG1SamplePoints));
    m_integratedNormalErrors = QVector<double>(numberOfCellEdges, 0);
    m_individualAngleErrors = QVector<QVector<double> >(numberOfCellEdges, QVector<double>(numberOfG1SamplePoints));
    m_integratedAngleErrors = QVector<double>(numberOfCellEdges, 0);

    m_individualBorderCurveLengths = QVector<double>(numberOfCellEdges);
    m_individualG1ErrorPositions = QVector<QVector<QVector3D> >(numberOfCellEdges, QVector<QVector3D>(numberOfG1SamplePoints));

    m_maxNormalError = -1;
    m_maxDetError = -1;
    m_maxAngleError = -1;

    m_totalIntegratedDetError = 0;
    m_totalIntegratedNormalError = 0;
    m_totalIntegratedAngleError = 0;
    m_totalBorderCurveLength = 0;

    //calculate G1 error
    #pragma omp parallel for
    for (int ceId = 0; ceId < numberOfCellEdges; ceId++) {
        CellEdge *ce = m_cellMesh->getEdge(ceId);
        const int incidentCellFaceIds[2] = {ce->getIncidentFaceId(0), ce->getIncidentFaceId(1)};
        CellFace *incidentCellFaces[2] = {m_cellMesh->getFace(incidentCellFaceIds[0]), m_cellMesh->getFace(incidentCellFaceIds[1])};
        const int internalCeId[2] = {incidentCellFaces[0]->getInternalEdgeId(ce), incidentCellFaces[1]->getInternalEdgeId(ce)};
        QVector<double> tangentLengths(numberOfG1SamplePoints);
        for (int iEval = 0; iEval < numberOfG1SamplePoints; iEval++) {
            double t = (double) iEval / (numberOfG1SamplePoints - 1);
            double u[2];
            double v[2];
            bool borderU[2];
            QVector3D Tu[2];
            QVector3D Tv[2];

            double length, backupLength;
            for (int k = 0; k < 2; k++) {
                double tOriented = t;
                if (incidentCellFaces[k]->edgeIsInverted(internalCeId[k]))
                    tOriented = 1 - t;
                if (internalCeId[k] == 0) {
                    u[k] = tOriented;
                    v[k] = 0;
                    borderU[k] = true;
                } else if (internalCeId[k] == 1) {
                    u[k] = 1;
                    v[k] = tOriented;
                    borderU[k] = false;
                } else if (internalCeId[k] == 2) {
                    u[k] = 1 - tOriented;
                    v[k] = 1;
                    borderU[k] = true;
                } else if (internalCeId[k] == 3) {
                    u[k] = 0;
                    v[k] = 1 - tOriented;
                    borderU[k] = false;
                }
                Tu[k] = m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceIds[k])->evaluateTu(u[k], v[k]);
                Tv[k] = m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceIds[k])->evaluateTv(u[k], v[k]);
                length = Tu[k].length();
                backupLength = Tv[k].length();
                Tu[k].normalize();
                Tv[k].normalize();
            }
            //u and v are still configured for the second surface
            m_individualG1ErrorPositions[ceId][iEval] = m_bSplinePatchNetwork->getBSplineSurface(incidentCellFaceIds[1])->evaluate(u[1], v[1]);

            QVector3D thirdVector(Tu[1]);
            if (borderU[1]) {
                thirdVector = Tv[1];
                length = backupLength;
            }
            tangentLengths[iEval] = length;

            Eigen::Matrix3d mat;
            mat << Tu[0].x(), Tv[0].x(), thirdVector.x(),
                   Tu[0].y(), Tv[0].y(), thirdVector.y(),
                   Tu[0].z(), Tv[0].z(), thirdVector.z();
            const double detErrorAbs = fabs(mat.determinant());
            QVector3D normals[2] = {QVector3D::crossProduct(Tu[0], Tv[0]).normalized(), QVector3D::crossProduct(Tu[1], Tv[1]).normalized()};
            const double normalErrorAbs = (normals[0] - normals[1]).length();
            double dotProduct = QVector3D::dotProduct(normals[0], normals[1]);
            if (dotProduct > 1)
                dotProduct = 1;
            else if (dotProduct < -1)
                dotProduct = -1;
            const double angleErrorAbs = fabs(acos(dotProduct) * 180/M_PI);

            m_individualDetErrors[ceId][iEval] = detErrorAbs;
            m_individualNormalErrors[ceId][iEval] = normalErrorAbs;
            m_individualAngleErrors[ceId][iEval] = angleErrorAbs;

            if (detErrorAbs > m_maxDetError)
                m_maxDetError = detErrorAbs;
            if (normalErrorAbs > m_maxNormalError)
                m_maxNormalError = normalErrorAbs;
            if (angleErrorAbs > m_maxAngleError)
                m_maxAngleError = angleErrorAbs;

            //integrate in parameter space
//            if (iEval > 0 && iEval < numberOfG1SamplePoints - 1) {
//                m_integratedDetErrors[ceId] += 2 * fabs(det);
//                m_integratedNormalErrors[ceId] += 2 * normalError;
//                m_integratedAngleErrors[ceId] += 2 * fabs(angleError);
//            } else {
//                m_integratedDetErrors[ceId] += fabs(det);
//                m_integratedNormalErrors[ceId] += normalError;
//                m_integratedAngleErrors[ceId] += fabs(angleError);
//            }
        }
//        m_integratedDetErrors[ceId] /= (2 * (numberOfG1SamplePoints - 1));
//        m_integratedNormalErrors[ceId] /= (2 * (numberOfG1SamplePoints - 1));
//        m_integratedAngleErrors[ceId] /= (2 * (numberOfG1SamplePoints - 1));

        //sum of errors times length of segments (i.e. neighborhoods)
//        double totalLength = 0;
//        QVector3D lastPoint = evalPositions[ceId][0];
//        QVector<double> segmentLengths(numberOfG1SamplePoints - 1);
//        for (int i = 1; i < numberOfG1SamplePoints; i++) {
//            const QVector3D currentPoint = evalPositions[ceId][i];
//            const double segmentLength = (lastPoint - currentPoint).length();
//            totalLength += segmentLength;
//            segmentLengths[i - 1] = segmentLength;
//            lastPoint = currentPoint;
//        }

//        double weightedTotalG1ErrorAngle = m_individualAngleErrors[ceId][0] * (segmentLengths[0]/2)
//                                         + m_individualAngleErrors[ceId][numberOfG1SamplePoints - 1] * (segmentLengths[numberOfG1SamplePoints - 2]/2);
//        double weightedTotalG1ErrorNormal = m_individualNormalErrors[ceId][0] * (segmentLengths[0]/2)
//                                          + m_individualNormalErrors[ceId][numberOfG1SamplePoints - 1] * (segmentLengths[numberOfG1SamplePoints - 2]/2);
//        double weightedTotalG1ErrorDet = m_individualDetErrors[ceId][0] * (segmentLengths[0]/2)
//                                       + m_individualDetErrors[ceId][numberOfG1SamplePoints - 1] * (segmentLengths[numberOfG1SamplePoints - 2]/2);

//        for (int i = 1; i < numberOfG1SamplePoints - 1; i++) {
//            const double factor = (segmentLengths[i - 1] + segmentLengths[i])/2;
//            weightedTotalG1ErrorAngle += m_individualAngleErrors[ceId][i] * factor;
//            weightedTotalG1ErrorNormal += m_individualNormalErrors[ceId][i] * factor;
//            weightedTotalG1ErrorDet += m_individualDetErrors[ceId][i] * factor;
//        }

        //integral of G1 error along curve in 3D space (trapezoidal rule)
        double integratedAngleError = 0;
        double integratedDetError = 0;
        double integratedNormalError = 0;
        for (int i = 0; i < numberOfG1SamplePoints - 1; i++) {
            integratedAngleError += (m_individualAngleErrors[ceId][i] * tangentLengths[i] + m_individualAngleErrors[ceId][i + 1] * tangentLengths[i + 1])/2;
            integratedDetError += (m_individualDetErrors[ceId][i] * tangentLengths[i] + m_individualDetErrors[ceId][i + 1] * tangentLengths[i + 1])/2;
            integratedNormalError += (m_individualNormalErrors[ceId][i] * tangentLengths[i] + m_individualNormalErrors[ceId][i + 1] * tangentLengths[i + 1])/2;
        }
        integratedAngleError /= (numberOfG1SamplePoints - 1);
        integratedDetError /= (numberOfG1SamplePoints - 1);
        integratedNormalError /= (numberOfG1SamplePoints - 1);

        m_integratedAngleErrors[ceId] = integratedAngleError;
        m_integratedNormalErrors[ceId] = integratedNormalError;
        m_integratedDetErrors[ceId] = integratedDetError;
        m_individualBorderCurveLengths[ceId] = m_bSplinePatchNetwork->getBSplineCurve(ceId)->calculateLengthNumerically(0, 1, numberOfG1SamplePoints);

    }

    for (int ceId = 0; ceId < numberOfCellEdges; ceId++) {
        m_totalIntegratedAngleError += m_integratedAngleErrors[ceId];
        m_totalIntegratedNormalError += m_integratedNormalErrors[ceId];
        m_totalIntegratedDetError += m_integratedDetErrors[ceId];
        m_totalBorderCurveLength += m_individualBorderCurveLengths[ceId];
    }

    m_totalAverageAngleError = m_totalIntegratedAngleError/m_totalBorderCurveLength;
    m_totalAverageNormalError = m_totalIntegratedNormalError/m_totalBorderCurveLength;
    m_totalAverageDetError = m_totalIntegratedDetError/m_totalBorderCurveLength;
}

void SplineEvaluator::evaluatePatchEnergyValues()
{
    const int numPatches = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
    m_patchEnergyValues.resize(numPatches);
    m_totalEnergy = 0;
    for (int i = 0; i < numPatches; i++) {
        const double energy = SplineEvaluator::calculatePatchEnergy(m_bSplinePatchNetwork->getBSplineSurface(i), *m_energyMatrix);
        m_patchEnergyValues[i] = energy;
        m_totalEnergy += energy;
    }
}

void SplineEvaluator::evaluatePatchAreas()
{
    m_patchAreasMesh = m_cellMesh->calculateApproximateSurfaceAreas();

    const int numberOfEvalPointsPerRow = 101;
    const int numPatches = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
    m_patchAreasBSplines.resize(numPatches);
    qDebug() << "TODO: test if surface area approximation with membrane energy is sufficient";

    m_totalAreaMesh = m_mesh->calculateTotalSurfaceArea();
    m_totalAreaBSplines = 0;

    #pragma omp parallel for
    for (int k = 0; k < numPatches; k++) {
        BSplineSurface *surface = m_bSplinePatchNetwork->getBSplineSurface(k);

        Vertex *v00 = new Vertex(QVector3D(0,0,0), -1);
        Vertex *v01 = new Vertex(QVector3D(0,0,0), -1);
        Vertex *v10 = new Vertex(QVector3D(0,0,0), -1);
        Vertex *v11 = new Vertex(QVector3D(0,0,0), -1);
        Edge *e0 = new Edge(v00, v10, -1);
        Edge *e1 = new Edge(v10, v11, -1);
        Edge *e2 = new Edge(v11, v01, -1);
        Edge *e3 = new Edge(v01, v00, -1);
        QVector<Edge *> edges;
        edges << e0 << e1 << e2 << e3;
        Face f(edges, QVector3D(0, 0, 1), -1);

        double area = 0;
        for (int i = 1; i < numberOfEvalPointsPerRow; i++) {
            for (int j = 1; j < numberOfEvalPointsPerRow; j++) {
                v00->updatePosition(surface->evaluate((double) (i-1) / (double) (numberOfEvalPointsPerRow-1), (double) (j-1) / (double) (numberOfEvalPointsPerRow-1)));
                v01->updatePosition(surface->evaluate((double) (i-1) / (double) (numberOfEvalPointsPerRow-1), (double) j / (double) (numberOfEvalPointsPerRow-1)));
                v10->updatePosition(surface->evaluate((double) i / (double) (numberOfEvalPointsPerRow-1), (double) (j-1) / (double) (numberOfEvalPointsPerRow-1)));
                v11->updatePosition(surface->evaluate((double) i / (double) (numberOfEvalPointsPerRow-1), (double) j / (double) (numberOfEvalPointsPerRow-1)));
                area += f.calculateAreaByCorners();
            }
        }

        m_patchAreasBSplines[k] = area;

        delete e0;
        delete e1;
        delete e2;
        delete e3;
        delete v00;
        delete v01;
        delete v10;
        delete v11;
    }

    for (int k = 0; k < numPatches; k++)
        m_totalAreaBSplines += m_patchAreasBSplines[k];
}

void SplineEvaluator::writeQualityFile(QString filename)
{
    const int numberOfControlPoints = m_bSplinePatchNetwork->getProperties()->getNumberOfControlPoints();
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();

    //Quality file (contains all relevant information)
    QFile fileQuality(filename);
    fileQuality.open(QIODevice::WriteOnly);
    QTextStream out(&fileQuality);



    out << "=== General Information ===\n"
        << "Number of border curves: " << QString::number(numberOfCellEdges) << "\n"
        << "Number of patches: " << QString::number(m_cellMesh->getNumberOfFaces()) << "\n"
        << "Number of data points: " << QString::number(m_mesh->getNumberOfVertices()) << "\n"
        << "Mesh size (X,Y,Z): " <<
           QString::number(m_meshSizeX) << " " <<
           QString::number(m_meshSizeY) << " " <<
           QString::number(m_meshSizeZ) << "\n"
        << "Bounding sphere (Cx, Cy, Cz): " <<
           QString::number(m_boundingSphereCenter.x()) << " " <<
           QString::number(m_boundingSphereCenter.y()) << " " <<
           QString::number(m_boundingSphereCenter.z()) << ", radius: " <<
           QString::number(m_boundingSphereRadius) << "\n"
        << "Average mesh edge length: " << QString::number(m_averageEdgeLength) << "\n"
        << "Number of control points per row: " << QString::number(numberOfControlPoints) << "\n"
        << "--- Properties ---\n"
        << "Area of the mesh: " << QString::number(m_totalAreaMesh)
        << " Area of the spline surface: " << QString::number(m_totalAreaBSplines) << "\n"
        << "Absolute Difference: " << QString::number(m_totalAreaMesh - m_totalAreaBSplines)
        << " Relative: " << QString::number(1 - m_totalAreaMesh/m_totalAreaBSplines) << "\n"
        << "Total energy: " << QString::number(m_totalEnergy) << " Energy/Area: " << QString::number(m_totalEnergy/m_totalAreaBSplines) << "\n"
        << "--- G1 error ---\n"
        << "Total determinant error: " << QString::number(m_totalIntegratedDetError) <<
           " Average: " << QString::number(m_totalAverageDetError) <<
           " max: " << QString::number(m_maxDetError) << "\n"
        << "Total normal error: " << QString::number(m_totalIntegratedNormalError) <<
           " Average: " << QString::number(m_totalAverageNormalError) <<
           " max: " << QString::number(m_maxNormalError) << "\n"
        << "Total angle error: " << QString::number(m_totalIntegratedAngleError) <<
           " Average: " << QString::number(m_totalAverageAngleError) <<
           " max: " << QString::number(m_maxAngleError) << "\n"
        << "Total length of border curves: " << QString::number(m_totalBorderCurveLength) << "\n"
        << "--- Fitting error ---\n"
        << "Total sum of squared Errors: " << QString::number(m_totalSquaredError) <<
           " rms error: " << QString::number(m_rmsError) <<
           " max error: " << QString::number(m_maxDataError) << "\n"
        << "rms/avgEdgeLength: " << QString::number(m_rmsError/m_averageEdgeLength) <<
           " maxError/avgEdgeLength: " << QString::number(m_maxDataError/m_averageEdgeLength) <<"\n"
        << "rms/maxMeshSize: " << QString::number(m_boundingBoxScaledRmsError) <<
           " maxError/maxSize: " << QString::number(m_boundingBoxScaledMaxDataError) <<"\n"
        << "RMS/boundingSphereDiameter [%]: " << QString::number(m_boundingSphereScaledRmsErrorPercent) <<
           " maxError/boundingSphereDiameter [%]: " << QString::number(m_boundingSphereScaledMaxDataErrorPercent) << "\n"
        << "===\n\n";

    for (int ceId = 0; ceId < numberOfCellEdges; ceId++) {
        double length = m_individualBorderCurveLengths[ceId];
        double errorDet = m_integratedDetErrors[ceId];
        double errorNorm = m_integratedNormalErrors[ceId];
        double errorAngle = m_integratedAngleErrors[ceId];
        out << "Border " << QString::number(ceId) << "\n";
        out << "Length: " << QString::number(length) << "\n";
        out << "Determinant error: " << QString::number(errorDet) << " Average: " << QString::number(errorDet/length) << "\n";
        out << "Normal error: " << QString::number(errorNorm) << " Average: " << QString::number(errorNorm/length) << "\n";
        out << "Angle error: " << QString::number(errorAngle) << " Average: " << QString::number(errorAngle/length) << "\n";
        out << "---\n";
    }
    fileQuality.close();
}

void SplineEvaluator::writeMatlabTriangulationStructureFile(QString filename)
{
    //each column contains the indices (i.e. line number) of the three vertices of one triangle
    QFile fileTri(filename);
    fileTri.open(QIODevice::WriteOnly);
    QTextStream outTri(&fileTri);
    const int numberOfTriangles = m_mesh->getNumberOfFaces();
    for (int i = 0; i < numberOfTriangles; i++) {
        Face *f = m_mesh->getFace(i);
        outTri << f->getVertex(0)->getId() + 1 << ", " << f->getVertex(1)->getId() + 1 << ", " << f->getVertex(2)->getId() + 1 << "\n"; //+1 because matlab starts counting at 1
    }
    fileTri.close();
}

void SplineEvaluator::writeMatlabRegularSamplingStructureFile(int samplePointsPerRow, QString filename)
{
    const int numberOfPatches = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
    const int samplePointsPerPatch = samplePointsPerRow * samplePointsPerRow;

    QFile fileStructure(filename);
    fileStructure.open(QIODevice::WriteOnly);
    QTextStream outStructure(&fileStructure);
    //Sample point indexing
    //
    // ...
    // ...
    // n+1 n+2 n+3
    // 1   2   3    ...  n

    for (int k = 0; k < numberOfPatches; k++) {
        for (int j = 0; j < samplePointsPerRow - 1; j++) {
            for (int i = 1; i < samplePointsPerRow; i++) {  //+1 because matlab starts counting at 1
                const int indexLowLeft = (k * samplePointsPerPatch) + (j * samplePointsPerRow) + i;
                const int indexLowRight = indexLowLeft + 1;
                const int indexUpLeft = indexLowLeft + samplePointsPerRow;
                const int indexUpRight = indexUpLeft + 1;

                outStructure << indexLowLeft << ", " << indexLowRight << ", " << indexUpRight << "\n"
                             << indexLowLeft << ", " << indexUpLeft << ", " << indexUpRight << "\n";
            }
        }
    }

    fileStructure.close();
}

void SplineEvaluator::writeMatlabRegularSampledCurvature(int samplePointsPerRow, QString filename)
{
    const int numberOfPatches = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();

    QFile fileCurvature(filename);
    fileCurvature.open(QIODevice::WriteOnly);
    QTextStream outCurvature(&fileCurvature);

    //Sample point indexing
    //
    // ...
    // ...
    // n+1 n+2 n+3
    // 1   2   3    ...  n
    for (int k = 0; k < numberOfPatches; k++) {
        BSplineSurface *surface = m_bSplinePatchNetwork->getBSplineSurface(k);
        for (int j = 0; j < samplePointsPerRow; j++) {
            for (int i = 0; i < samplePointsPerRow; i++) {
                const double u = (double) i / (samplePointsPerRow - 1);
                const double v = (double) j / (samplePointsPerRow - 1);

                double K, H;    //K - gaussian, H - mean curvature
                surface->evaluateGaussianAndMeanCurvature(u, v, K, H);
                const double k1k1_k2k2 = ((4 * H * H) - (2 * K));  //equivalent to k1*k1 + k2*k2, since H = (k1+k2)/2, K = k1*k2
                const QVector3D surfUU = surface->evaluateDerivative(u, v, 2, 0);
                const QVector3D surfUV = surface->evaluateDerivative(u, v, 1, 1);
                const QVector3D surfVV = surface->evaluateDerivative(u, v, 0, 2);
                const double thinPlate = surfUU.lengthSquared() + 2 * surfUV.lengthSquared() + surfVV.lengthSquared();

                const QVector3D bSplinePos = surface->evaluate(u, v);

                outCurvature << bSplinePos.x() << ", " << bSplinePos.y() << ", " << bSplinePos.z() << ", " << K << ", " << H << ", " << k1k1_k2k2 << ", " << thinPlate << "\n";
            }
        }
    }

    fileCurvature.close();
}

void SplineEvaluator::writeMatlabRegularSamplingG1Error(QString filename)
{
    QFile g1File(filename);
    g1File.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&g1File);

    const int numCellEdges = m_bSplinePatchNetwork->getNumberOfBSplineCurves();
    const int numSamplePoints = m_individualAngleErrors[0].size();

    //first line: numberOfCellEdges, numberOfSamplePointsPerEdge
    out << numCellEdges << ", " << numSamplePoints << ", -1, -1, -1, -1\n";   //fill with dummy values so that matlab can load it as matrix
    //posX, posY, posZ, errorNormal, errorAngle
    for (int iCe = 0; iCe < numCellEdges; iCe++) {
        for (int i = 0; i < numSamplePoints; i++) {
            QVector3D pos = m_individualG1ErrorPositions[iCe][i];
            out << pos.x() << ", " << pos.y() << ", " << pos.z() << ", "
                << m_individualNormalErrors[iCe][i] << ", "
                << m_individualAngleErrors[iCe][i] << "\n";
        }
    }

    g1File.close();
}

void SplineEvaluator::writeMatlabDataErrorFiles(QString filenameTri, QString filenameBSpline)
{
    const int numberOfDataPoints = m_mesh->getNumberOfVertices();

    //column 1,2,3: triangle vertices, column 4: data error
    QFile fileTriError(filenameTri);
    //column 1,2,3: projection of triangle vertices onto B-Spline surface, column 4: data error
    QFile fileBSplineError(filenameBSpline);
    fileTriError.open(QIODevice::WriteOnly);
    fileBSplineError.open(QIODevice::WriteOnly);
    QTextStream outTriError(&fileTriError);
    QTextStream outBSplineError(&fileBSplineError);
    for (int i = 0; i < numberOfDataPoints; i++) {
        Vertex *p = m_mesh->getVertex(i);
        outTriError << p->getPosX() << ", " << p->getPosY() << ", " << p->getPosZ() << ", " << m_dataErrors[i] << ", " << m_dataErrors[i]/(m_boundingSphereRadius * 2) * 100 << "\n";

        const QVector3D bSplinePos = m_projectedCoordinates[i];
        outBSplineError << bSplinePos.x() << ", " << bSplinePos.y() << ", " << bSplinePos.z() << ", " << m_dataErrors[i] << ", " << m_dataErrors[i]/(m_boundingSphereRadius * 2) * 100 << "\n";

    }
    fileTriError.close();
    fileBSplineError.close();
}

void SplineEvaluator::writeMatlabG1ErrorAtPLVFiles(QString filenameTemplateTri, QString filenameTemplateBSpline)
{
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();

    //column 1,2,3: polyline vertex, column 4: G1 error
    //column 1,2,3: projection of plv onto B-Spline surface, column 4, 5: G1 error (normal difference, angles)
    for (int iCE = 0; iCE < numberOfCellEdges; iCE++) {
        QFile fileTriG1Error(filenameTemplateTri.arg(iCE));
        QFile fileBSG1Error(filenameTemplateBSpline.arg(iCE));
        fileTriG1Error.open(QIODevice::WriteOnly);
        fileBSG1Error.open(QIODevice::WriteOnly);
        QTextStream outTriG1Error(&fileTriG1Error);
        QTextStream outBSG1Error(&fileBSG1Error);

        CellEdge *ce = m_cellMesh->getEdge(iCE);
        const int numPLV = ce->getPolylineVertices()->size();
        const int fId0 = ce->getIncidentFaceId(0);
        const int fId1 = ce->getIncidentFaceId(1);
        Parameterization *param0 = m_cellMesh->getQuadCellParameterization(fId0);
        Parameterization *param1 = m_cellMesh->getQuadCellParameterization(fId1);
        for (int i = 0; i < numPLV; i++) {
            PolylineVertex *plv = ce->getPolylineVertices()->at(i);
            const QVector2D paramValue0 = param0->getParameter(plv);
            double errorNormal = -1;
            double errorAngle = -1;
            if (ce->isG1Edge()) {
                const QVector2D paramValue1 = param1->getParameter(plv);
                const QVector3D Tu0 = m_bSplinePatchNetwork->getBSplineSurface(fId0)->evaluateTu(paramValue0.x(), paramValue0.y());
                const QVector3D Tv0 = m_bSplinePatchNetwork->getBSplineSurface(fId0)->evaluateTv(paramValue0.x(), paramValue0.y());
                const QVector3D Tu1 = m_bSplinePatchNetwork->getBSplineSurface(fId1)->evaluateTu(paramValue1.x(), paramValue1.y());
                const QVector3D Tv1 = m_bSplinePatchNetwork->getBSplineSurface(fId1)->evaluateTv(paramValue1.x(), paramValue1.y());
                const QVector3D normal0 = QVector3D::crossProduct(Tu0, Tv0).normalized();
                const QVector3D normal1 = QVector3D::crossProduct(Tu1, Tv1).normalized();
                errorNormal = (normal0 - normal1).length();
                errorAngle = acos((2 - errorNormal * errorNormal)/2) * 180/M_PI;
            }
            outTriG1Error << plv->getPosX() << ", " << plv->getPosY() << ", " << plv->getPosZ() << ", " << errorNormal << ", " << errorAngle << "\n";

            const QVector3D bsplinePos = m_bSplinePatchNetwork->getBSplineSurface(fId0)->evaluate(paramValue0.x(), paramValue0.y());
            outBSG1Error << bsplinePos.x() << ", " << bsplinePos.y() << ", " << bsplinePos.z() << ", " << errorNormal << ", " << errorAngle << "\n";
        }

        fileTriG1Error.close();
        fileBSG1Error.close();
    }
}

double SplineEvaluator::getRmsError()
{
    return m_rmsError;
}

double SplineEvaluator::getMaxDataError()
{
    return m_maxDataError;
}

double SplineEvaluator::getTotalAverageDetError()
{
    return m_totalAverageDetError;
}

double SplineEvaluator::getTotalAverageNormalError()
{
    return m_totalAverageNormalError;
}

double SplineEvaluator::getTotalAverageAngleError()
{
    return m_totalAverageAngleError;
}

double SplineEvaluator::getMaxDetError()
{
    return m_maxDetError;
}

double SplineEvaluator::getMaxNormalError()
{
    return m_maxNormalError;
}

double SplineEvaluator::getMaxAngleError()
{
    return m_maxAngleError;
}

double SplineEvaluator::getBoundingBoxScaledRmsError()
{
    return m_boundingBoxScaledRmsError;
}

double SplineEvaluator::getBoundingBoxScaledMaxDataError()
{
    return m_boundingBoxScaledMaxDataError;
}

double SplineEvaluator::getBoundingSphereScaledRmsError()
{
    return m_boundingSphereScaledRmsErrorPercent;
}

double SplineEvaluator::getBoundingSphereScaledMaxDataError()
{
    return m_boundingSphereScaledMaxDataErrorPercent;
}

double SplineEvaluator::getTotalSurfaceAreaMesh()
{
    return m_totalAreaMesh;
}

double SplineEvaluator::getTotalSurfaceAreaBSplines()
{
    return m_totalAreaBSplines;
}

double SplineEvaluator::getTotalEnergy()
{
    return m_totalEnergy;
}

double SplineEvaluator::calculatePatchEnergy(BSplineSurface *surface, Eigen::SparseMatrix<double> &energyMatrix)
{
    const int numberOfControlPointsPerRow = surface->getControlPointsP()->size();
    const int numberOfControlPointsPerPatch = numberOfControlPointsPerRow * numberOfControlPointsPerRow;
    Eigen::VectorXd ctrX(numberOfControlPointsPerPatch);
    Eigen::VectorXd ctrY(numberOfControlPointsPerPatch);
    Eigen::VectorXd ctrZ(numberOfControlPointsPerPatch);
    QVector<QVector3D> ctrPointsAsVector = surface->getControlPointsAsVector();

    for (int i = 0; i < numberOfControlPointsPerPatch; i++) {
        ctrX[i] = ctrPointsAsVector[i].x();
        ctrY[i] = ctrPointsAsVector[i].y();
        ctrZ[i] = ctrPointsAsVector[i].z();
    }

    const double energyX = ctrX.transpose() * energyMatrix * ctrX;
    const double energyY = ctrY.transpose() * energyMatrix * ctrY;
    const double energyZ = ctrZ.transpose() * energyMatrix * ctrZ;

    return energyX + energyY + energyZ;
}

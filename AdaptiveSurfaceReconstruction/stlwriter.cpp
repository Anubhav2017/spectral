#include "stlwriter.h"

#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QDir>

#include "stdint.h"
#include <iostream>
#include <fstream>

#include <QDebug>

StlWriter::StlWriter(QString defaultFolder, bool cleanDirectoryBeforeWriting, bool askOverwrite, bool consoleOutput):
    m_defaultFolder(defaultFolder), m_askOverwrite(askOverwrite), m_consoleOutput(consoleOutput)
{
    //QDir(defaultFolder).removeRecursively();  //Only works with Qt 5 or newer
    if (cleanDirectoryBeforeWriting)
        removeRecursively(defaultFolder);

    if (!QDir(defaultFolder).exists())
        QDir().mkpath(defaultFolder);

    QDir().mkpath(defaultFolder+"/debug/");
}

void StlWriter::copyAllFilesTo(QString destinationFolder)
{
    bool worked = copyRecursively(m_defaultFolder, destinationFolder);
    if (m_consoleOutput && worked)
        std::cout << QString("Files copied to %1").arg(destinationFolder).toStdString();
    else if (!worked)
        std::cout << QString("Error: Could not copy files from %1 to %2").arg(m_defaultFolder).arg(destinationFolder).toStdString();

}

void StlWriter::writeMeshFacesToStl(Mesh *m, QString filename)
{
    const int numFaces = m->getNumberOfFaces();
    QVector<Face *> faces(numFaces);

    for (int i = 0; i < numFaces; i++)
        faces[i] = m->getFace(i);

    //this->writeFacesToStl(&faces, filename);
    this->writeFacesToStlBinary(&faces, filename);
}

#include "stdint.h"
void StlWriter::writeFacesToStlBinary(QVector<Face *> *faces, QString filename)
{
    QString fullFilename = QString(m_defaultFolder).append(filename);
    if (!fullFilename.endsWith(".stl", Qt::CaseInsensitive))
        fullFilename.append(".stl");

    std::ofstream stream;
    std::string stdFilename = fullFilename.toStdString();
    const int n = stdFilename.size();
    char name[n + 1];
    std::strcpy(name, stdFilename.c_str());
    stream.open(name, std::ios::binary | std::ios::out);

    std::string headstring = "Binary STL File";
    char head[80] = {0};
    strcpy(head, headstring.c_str());
    stream.write(head, 80);

    //TODO find types that work under windows and linux
    const uint32_t numberOfFaces = faces->size();
    stream.write((char*) &numberOfFaces, 4);

    //char attribute[2]="0";
    uint16_t attCount = 0;

    for (int i = 0; i < faces->size(); i++) {
        Face *f = faces->at(i);

        const QVector3D v0 = f->getVertex(0)->getPosition();
        float corner0[3] = {(float) v0.x(), (float) v0.y(), (float) v0.z()};
        float normal[3] = {(float) f->getNormal().x(), (float) f->getNormal().y(), (float) f->getNormal().z()};

        for (int j = 2; j < f->getNumberOfVertices(); j++) {
            QVector3D v1 = f->getVertex(j - 1)->getPosition();
            QVector3D v2 = f->getVertex(j)->getPosition();

            float corner1[3] = {(float) v1.x(), (float) v1.y(), (float) v1.z()};
            float corner2[3] = {(float) v2.x(), (float) v2.y(), (float) v2.z()};

            stream.write((char *) &normal[0], 4);
            stream.write((char *) &normal[1], 4);
            stream.write((char *) &normal[2], 4);
            stream.write((char *) &corner0[0], 4);
            stream.write((char *) &corner0[1], 4);
            stream.write((char *) &corner0[2], 4);
            stream.write((char *) &corner1[0], 4);
            stream.write((char *) &corner1[1], 4);
            stream.write((char *) &corner1[2], 4);
            stream.write((char *) &corner2[0], 4);
            stream.write((char *) &corner2[1], 4);
            stream.write((char *) &corner2[2], 4);

            stream.write((char *) &attCount, 2);
        }
    }

    stream.close();
}

void StlWriter::writeFacesToStlASCII(QVector<Face *> *faces, QString filename)
{
    QFile *file = processFilename(filename);
    if (!file)
        return;

    file->open(QIODevice::WriteOnly);
    QTextStream out(file);

    out << "solid Triangulation\n";
    for (int i = 0; i < faces->size(); i++) {
        Face *f = faces->at(i);

        QVector3D v0 = f->getVertex(0)->getPosition();
        QVector3D normal = f->getNormal();

        for (int j = 2; j < f->getNumberOfVertices(); j++) {
            QVector3D v1 = f->getVertex(j - 1)->getPosition();
            QVector3D v2 = f->getVertex(j)->getPosition();

            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(v0.x()) << " "
                << QString::number(v0.y()) << " "
                << QString::number(v0.z()) << "\n"
                << "      vertex "
                << QString::number(v1.x()) << " "
                << QString::number(v1.y()) << " "
                << QString::number(v1.z()) << "\n"
                << "      vertex "
                << QString::number(v2.x()) << " "
                << QString::number(v2.y()) << " "
                << QString::number(v2.z()) << "\n"
                << "    endloop\n  endfacet\n";
        }
    }

    out << "endsolid Triangulation";

    file->close();
    delete file;
    if (m_consoleOutput)
        std::cout << "done." << std::endl;
}

void StlWriter::writeEdgesToStl(QVector<Edge *> *edges, QString filename, double width)
{
    QFile *file = processFilename(filename);
    if (!file)
        return;

    file->open(QIODevice::WriteOnly);
    QTextStream out(file);

    out << "solid Triangulation\n";
    for (int i = 0; i < edges->size(); i++) {
        Edge *e = edges->at(i);
        writeLineSegment(e->getVertex(0)->getPosition(), e->getVertex(1)->getPosition(), &out, width);
    }

    out << "endsolid Triangulation";

    file->close();
    delete file;
    if (m_consoleOutput)
        std::cout << "done." << std::endl;
}

void StlWriter::writeBSplineCurvesToStl(QVector<BSplineCurve *> curves, QString filename, double linewidth, int numberOfEvalPointsPerCurve)
{
    const int numberOfCurves = curves.size();

    QVector<Edge *> tmpEdges;
    for (int iCurve = 0; iCurve < numberOfCurves; iCurve++) {
        BSplineCurve *curve = curves[iCurve];
        for (int i = 1; i < numberOfEvalPointsPerCurve; i++) {
            Vertex *v0 = new Vertex(curve->evaluateDeBoor((double) (i - 1) / (double) (numberOfEvalPointsPerCurve-1)), -1);
            Vertex *v1 = new Vertex(curve->evaluateDeBoor((double)  i      / (double) (numberOfEvalPointsPerCurve-1)), -1);
            tmpEdges.append(new Edge(v0, v1, -1));
        }
    }
    this->writeEdgesToStl(&tmpEdges, filename, linewidth);

    //clean up
    for (int i = 0; i < tmpEdges.size(); i++) {
        delete tmpEdges[i]->getVertex(0);
        delete tmpEdges[i]->getVertex(1);
        delete tmpEdges[i];
    }
}

void StlWriter::writeBSplineSurfacesToStl(QVector<BSplineSurface *> surfaces, QString filename, int numberOfEvalPointsPerRow)
{
    const int numberOfSurfaces = surfaces.size();
    QVector<Face *> tmpFaces;
    for (int iSurf = 0; iSurf < numberOfSurfaces; iSurf++) {
        BSplineSurface *surface = surfaces[iSurf];
        for (int i = 1; i < numberOfEvalPointsPerRow; i++) {
            for (int j = 1; j < numberOfEvalPointsPerRow; j++) {
                Vertex *v00 = new Vertex(surface->evaluate((double) (i-1) / (double) (numberOfEvalPointsPerRow-1), (double) (j-1) / (double) (numberOfEvalPointsPerRow-1)), -1);
                Vertex *v01 = new Vertex(surface->evaluate((double) (i-1) / (double) (numberOfEvalPointsPerRow-1), (double) j / (double) (numberOfEvalPointsPerRow-1)), -1);
                Vertex *v10 = new Vertex(surface->evaluate((double) i / (double) (numberOfEvalPointsPerRow-1), (double) (j-1) / (double) (numberOfEvalPointsPerRow-1)), -1);
                Vertex *v11 = new Vertex(surface->evaluate((double) i / (double) (numberOfEvalPointsPerRow-1), (double) j / (double) (numberOfEvalPointsPerRow-1)), -1);
                Edge *e0 = new Edge(v00, v10, -1);
                Edge *e1 = new Edge(v10, v11, -1);
                Edge *e2 = new Edge(v11, v01, -1);
                Edge *e3 = new Edge(v01, v00, -1);
                QVector<Edge *> edges;
                edges << e0 << e1 << e2 << e3;
                tmpFaces << new Face(edges, QVector3D(0, 0, 1), -1);
            }
        }
    }

    this->writeFacesToStlASCII(&tmpFaces, filename);

    //clean up
    const int numTmpFaces = tmpFaces.size();
    for (int i = 0; i < numTmpFaces; i++) {
        for (int j = 0; j < 4; j++) {
            delete tmpFaces[i]->getVertex(j);
            delete tmpFaces[i]->getEdge(j);
        }
        delete tmpFaces[i];
    }
}

void StlWriter::writeBSplineSubdivisionsToStl(QVector<BSplineSurface *> surfaces, QVector<SubSurfaceTree *> subdivisions, QString filename, double linewidth, int numberOfEvalPoints)
{
    const int numSurfaces = surfaces.size();

    if (subdivisions.size() != numSurfaces)
        return;

    QVector<Vertex *> tmpVertices;
    QVector<Edge *> tmpEdges;
    for (int iSurf = 0; iSurf < numSurfaces; iSurf++) {
        BSplineSurface *surf = surfaces[iSurf];
        SubSurfaceTree *root = subdivisions[iSurf];
        QVector<SubSurfaceTree *> leaves = root->getAllLeaves();
        const int numLeaves = leaves.size();
        const double uMax = root->getU1();
        const double vMax = root->getV1();

        for (int iLeaf = 0; iLeaf < numLeaves; iLeaf++) {
            SubSurfaceTree *leaf = leaves[iLeaf];
            const double u0 = leaf->getU0();
            const double u1 = leaf->getU1();
            const double v0 = leaf->getV0();
            const double v1 = leaf->getV1();
            const double uDiff = u1 - u0;
            const double vDiff = v1 - v0;

            //line (u0, v1) to (u1, v1)
            if (v1 != vMax) {
                Vertex *vLast = new Vertex (surf->evaluate(u0, v1), -1);
                tmpVertices << vLast;
                for (int i = 1; i < numberOfEvalPoints; i++) {
                    const double u = u0 + uDiff * (double) i / (numberOfEvalPoints - 1);
                    Vertex *vNext = new Vertex(surf->evaluate(u, v1), -1);
                    Edge *e = new Edge(vLast, vNext, -1);
                    tmpEdges << e;
                    tmpVertices << vNext;
                    vLast = vNext;
                }
            }

            //line (u1, v0) to (u1, v1)
            if (u1 != uMax) {
                Vertex *vLast = new Vertex (surf->evaluate(u1, v0), -1);
                tmpVertices << vLast;
                for (int i = 1; i < numberOfEvalPoints; i++) {
                    const double v = v0 + vDiff * (double) i / (numberOfEvalPoints - 1);
                    Vertex *vNext = new Vertex(surf->evaluate(u1, v), -1);
                    Edge *e = new Edge(vLast, vNext, -1);
                    tmpEdges << e;
                    tmpVertices << vNext;
                    vLast = vNext;
                }
            }
        }
    }

    this->writeEdgesToStl(&tmpEdges, filename, linewidth);

    const int numTmpEdges = tmpEdges.size();
    for (int i = 0; i < numTmpEdges; i++)
        delete tmpEdges[i];

    const int numTmpVertices = tmpVertices.size();
    for (int i = 0; i < numTmpVertices; i++)
        delete tmpVertices[i];
}

void StlWriter::writeCellFaceContentToStl(CellFace *cf, QString filename, double vertexSize, double edgeWidth)
{
    QFile *file = processFilename(filename);
    if (!file)
        return;

    QVector<QVector3D> points;
    foreach (Vertex *v, *cf->getMeshVertices()) {
        points.append(v->getPosition());
    }

    file->open(QIODevice::WriteOnly);
    QTextStream out(file);

    out << "solid Triangulation\n";
    for (int i = 0; i < points.size(); i++) {
        QVector3D vCenter = points.at(i);

        QVector3D v[4];
        v[0] = vCenter + QVector3D(vertexSize, 0, 0);
        v[1] = vCenter + QVector3D(0, vertexSize, 0);
        v[2] = vCenter + QVector3D(-vertexSize, 0, 0);
        v[3] = vCenter + QVector3D(0, -vertexSize, 0);
        QVector3D vTop = vCenter + QVector3D(0, 0, vertexSize);
        QVector3D vBottom = vCenter + QVector3D(0, 0, -vertexSize);

        for (int j = 0; j < 4; j++) {
            QVector3D vCurrent = v[j];
            QVector3D vNext = v[(j+1)%4];
            QVector3D normal = QVector3D::crossProduct(vNext-vCurrent, vTop-vCurrent);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "      vertex "
                << QString::number(vTop.x()) << " "
                << QString::number(vTop.y()) << " "
                << QString::number(vTop.z()) << "\n"
                << "    endloop\n  endfacet\n";

            normal = QVector3D::crossProduct(vNext-vCurrent, vCurrent-vBottom);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vBottom.x()) << " "
                << QString::number(vBottom.y()) << " "
                << QString::number(vBottom.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "    endloop\n  endfacet\n";
        }
    }

    for (int i = 0; i < cf->getMeshEdges()->size(); i++) {
        Edge *e = cf->getMeshEdges()->at(i);
        writeLineSegment(e->getVertex(0)->getPosition(), e->getVertex(1)->getPosition(), &out, edgeWidth);
    }

    out << "endsolid Triangulation";

    file->close();
    delete file;
    if (m_consoleOutput)
        std::cout << "done." << std::endl;
}

void StlWriter::writeCellBorderContentToStl(CellFace *cf, Mesh *originalMesh, QString filename, double vertexSize, double edgeWidth)
{
    QFile *file = processFilename(filename);
    if (!file)
        return;

    file->open(QIODevice::WriteOnly);
    QTextStream out(file);

    out << "solid Triangulation\n";

    QVector<QVector3D> points;
    for (int i = 0; i < cf->getNumberOfEdges(); i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices())
            if (plv->isMeshVertex()) {
                points.append(plv->getMeshVertex()->getPosition());
                foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
                    Edge *e = originalMesh->getEdge(eId);
                    writeLineSegment(e->getVertex(0)->getPosition(), e->getVertex(1)->getPosition(), &out, edgeWidth);
                }
            } else if (plv->isCrossingVertex()) {
                Edge *e = plv->getCrossingEdge();
                writeLineSegment(e->getVertex(0)->getPosition(), e->getVertex(1)->getPosition(), &out, edgeWidth);
            }
    }

    for (int i = 0; i < points.size(); i++) {
        QVector3D vCenter = points.at(i);

        QVector3D v[4];
        v[0] = vCenter + QVector3D(vertexSize, 0, 0);
        v[1] = vCenter + QVector3D(0, vertexSize, 0);
        v[2] = vCenter + QVector3D(-vertexSize, 0, 0);
        v[3] = vCenter + QVector3D(0, -vertexSize, 0);
        QVector3D vTop = vCenter + QVector3D(0, 0, vertexSize);
        QVector3D vBottom = vCenter + QVector3D(0, 0, -vertexSize);

        for (int j = 0; j < 4; j++) {
            QVector3D vCurrent = v[j];
            QVector3D vNext = v[(j+1)%4];
            QVector3D normal = QVector3D::crossProduct(vNext-vCurrent, vTop-vCurrent);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "      vertex "
                << QString::number(vTop.x()) << " "
                << QString::number(vTop.y()) << " "
                << QString::number(vTop.z()) << "\n"
                << "    endloop\n  endfacet\n";

            normal = QVector3D::crossProduct(vNext-vCurrent, vCurrent-vBottom);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vBottom.x()) << " "
                << QString::number(vBottom.y()) << " "
                << QString::number(vBottom.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "    endloop\n  endfacet\n";
        }
    }

    out << "endsolid Triangulation";

    file->close();
    delete file;
    if (m_consoleOutput)
        std::cout << "done." << std::endl;
}

void StlWriter::writeControlPointsOfIndividualPhases(BSplinePatchNetwork *patchNetwork, QString filenameTemplate, double vertexSize)
{
    QVector<QVector3D> cornerPoints, edgePoints, interiorPoints;

    const int numberOfPatches = patchNetwork->getNumberOfBSplineSurfaces();
    const int numberOfCtrPoints = patchNetwork->getProperties()->getNumberOfControlPoints();

    for (int iPatch = 0; iPatch < numberOfPatches; iPatch++) {
        QVector<QVector<QVector3D> > *ctrPoints = patchNetwork->getBSplineSurface(iPatch)->getControlPointsP();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                cornerPoints << ctrPoints->at(i).at(j);
                cornerPoints << ctrPoints->at(numberOfCtrPoints - 1 - i).at(j);
                cornerPoints << ctrPoints->at(i).at(numberOfCtrPoints - 1 - j);
                cornerPoints << ctrPoints->at(numberOfCtrPoints - 1 - i).at(numberOfCtrPoints - 1 - j);
            }
        }

        for (int i = 2; i < numberOfCtrPoints - 2; i++) {
            for (int j = 0; j < 2; j++) {
                edgePoints << ctrPoints->at(i).at(j);
                edgePoints << ctrPoints->at(i).at(numberOfCtrPoints - 1 - j);
                edgePoints << ctrPoints->at(j).at(i);
                edgePoints << ctrPoints->at(numberOfCtrPoints - 1 - j).at(i);
            }
        }

        for (int i = 2; i < numberOfCtrPoints - 2; i++)
            for(int j = 2; j < numberOfCtrPoints - 2; j++)
                interiorPoints << ctrPoints->at(i).at(j);
    }

    writePoints(&cornerPoints, vertexSize, filenameTemplate.arg("corners"));
    writePoints(&edgePoints, vertexSize, filenameTemplate.arg("edges"));
    writePoints(&interiorPoints, vertexSize, filenameTemplate.arg("interior"));
}

void StlWriter::writeBSplineControlPointsToStl(BSplinePatchNetwork *patchNetwork, QString filenameEdges, QString filenameVertices, double vertexSize, double linewidth)
{
    QVector<Edge *> controlPolygonEdges;
    QVector<QVector3D> controlPointVertices;
    const int numberOfSurfaces = patchNetwork->getNumberOfBSplineSurfaces();
    for (int k = 0; k < numberOfSurfaces; k++) {
        BSplineSurface *surface = patchNetwork->getBSplineSurface(k);
        const int numberOfControlPointsU = surface->getPropertiesU()->getNumberOfControlPoints();
        const int numberOfControlPointsV = surface->getPropertiesU()->getNumberOfControlPoints();
        QVector<QVector<QVector3D> > *controlPoints = surface->getControlPointsP();
        QVector<QVector3D> vertices;
        QVector<Edge *> edges;
        for (int i = 0; i < numberOfControlPointsU; i++) {
            vertices << controlPoints->at(i);
            for (int j = 0; j < numberOfControlPointsV; j++) {
                if (i > 0) {
                    Vertex *v0 = new Vertex(controlPoints->at(i-1).at(j), -1);
                    Vertex *v1 = new Vertex(controlPoints->at(i).at(j), -1);
                    Edge *e = new Edge(v0, v1, -1);
                    edges << e;
                }
                if (j > 0) {
                    Vertex *v0 = new Vertex(controlPoints->at(i).at(j-1), -1);
                    Vertex *v1 = new Vertex(controlPoints->at(i).at(j), -1);
                    Edge *e = new Edge(v0, v1, -1);
                    edges << e;
                }
            }
        }
        controlPolygonEdges << edges;
        controlPointVertices << vertices;
    }

    this->writeEdgesToStl(&controlPolygonEdges, filenameEdges, linewidth);
    this->writePoints(&controlPointVertices, vertexSize, filenameVertices);

    //clean up
    for (int i = 0; i < controlPolygonEdges.size(); i++) {
        delete controlPolygonEdges[i]->getVertex(0);
        delete controlPolygonEdges[i]->getVertex(1);
        delete controlPolygonEdges[i];
    }
}

void StlWriter::writeBSplineIsolinesToStl(QVector<BSplineSurface *> surfaces, QVector<double> isolinesU, QVector<double> isolinesV, QString filename, QPair<double, double> parameterInterval, int numberOfEvalPoints, double linewidth)
{
    QVector<Vertex *> vertices;
    QVector<Edge *> edges;
    foreach (BSplineSurface *surf, surfaces) {
        foreach (double u, isolinesU) {
            const int offset = vertices.size();
            for (int i = 0; i < numberOfEvalPoints; i++) {
                const double t = parameterInterval.first + parameterInterval.second * (double) i / (numberOfEvalPoints - 1);
                vertices << new Vertex(surf->evaluate(u, t), -1);
            }

            for (int i = 1; i < numberOfEvalPoints; i++)
                edges << new Edge(vertices[offset + i - 1], vertices[offset + i], -1);
        }

        foreach (double v, isolinesV) {
            const int offset = vertices.size();
            for (int i = 0; i < numberOfEvalPoints; i++) {
                const double t = parameterInterval.first + parameterInterval.second * (double) i / (numberOfEvalPoints - 1);
                vertices << new Vertex(surf->evaluate(t, v), -1);
            }

            for (int i = 1; i < numberOfEvalPoints; i++)
                edges << new Edge(vertices[offset + i - 1], vertices[offset + i], -1);
        }
    }

    this->writeEdgesToStl(&edges, filename, linewidth);

    //clean up
    const int numVertices = vertices.size();
    const int numEdges = edges.size();

    for (int i = 0; i < numEdges; i++)
        delete edges[i];

    for (int i = 0; i < numVertices; i++)
        delete vertices[i];
}

void StlWriter::writeBSplineDataErrorToStl(BSplinePatchNetwork *patchNetwork, QString filename, QString indSurfaceNameTemplate, bool errorOnBoundaries, double edgeWidth)
{
    const int numberOfSurfaces = patchNetwork->getNumberOfBSplineSurfaces();
    bool doIndividualDataErrors = !indSurfaceNameTemplate.isEmpty();
    CellMesh *cm = patchNetwork->getCellMesh();

    QVector<Edge*> allDataErrors;
    for (int cfId = 0; cfId < numberOfSurfaces; cfId++) {
        QVector<Edge*> dataErrors;
        CellFace *cf = cm->getFace(cfId);
        Parameterization *cfParam = cm->getQuadCellParameterization(cfId);
        BSplineSurface *surf = patchNetwork->getBSplineSurface(cfId);
        const int numberOfDataPoints = cf->getMeshVertices()->size();
        for (int i = 0; i < numberOfDataPoints; i++) {
            Vertex *data = cf->getMeshVertices()->at(i);
            QVector2D param = cfParam->getParameter(data);
            Vertex *interpolation = new Vertex(surf->evaluate(param.x(), param.y()), -1);
            dataErrors << new Edge(data, interpolation, -1);
        }

        if (errorOnBoundaries) {
            for (int eId = 0; eId < 4; eId++) {
                CellEdge *ce = (CellEdge *) cf->getEdge(eId);
                const int numPLV = ce->getPolylineVertices()->size();
                for (int i = 0; i < numPLV; i++) {
                    if (ce->getPolylineVertices()->at(i)->isMeshVertex()) {
                        Vertex *data = ce->getPolylineVertices()->at(i);
                        QVector2D param = cfParam->getParameter(ce->getPolylineVertices()->at(i));
                        Vertex *interpolation = new Vertex(surf->evaluate(param.x(), param.y()), -1);
                        dataErrors << new Edge(data, interpolation, -1);
                    }
                }
            }
        }

        if (doIndividualDataErrors)
            this->writeEdgesToStl(&dataErrors, indSurfaceNameTemplate.arg(cfId), edgeWidth);

        allDataErrors << dataErrors;
    }
    this->writeEdgesToStl(&allDataErrors, filename, edgeWidth);

    //clean up
    for (int i = 0; i < allDataErrors.size(); i++) {
        delete allDataErrors[i]->getVertex(1);
        delete allDataErrors[i];
    }
}

void StlWriter::writePoints(QVector<QVector3D> *points, double size, QString filename)
{
    QFile *file = processFilename(filename);
    if (!file)
        return;

    file->open(QIODevice::WriteOnly);
    QTextStream out(file);

    out << "solid Triangulation\n";
    for (int i = 0; i < points->size(); i++) {
        QVector3D vCenter = points->at(i);

        QVector3D v[4];
        v[0] = vCenter + QVector3D(size, 0, 0);
        v[1] = vCenter + QVector3D(0, size, 0);
        v[2] = vCenter + QVector3D(-size, 0, 0);
        v[3] = vCenter + QVector3D(0, -size, 0);
        QVector3D vTop = vCenter + QVector3D(0, 0, size);
        QVector3D vBottom = vCenter + QVector3D(0, 0, -size);

        for (int j = 0; j < 4; j++) {
            QVector3D vCurrent = v[j];
            QVector3D vNext = v[(j+1)%4];
            QVector3D normal = QVector3D::crossProduct(vNext-vCurrent, vTop-vCurrent);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "      vertex "
                << QString::number(vTop.x()) << " "
                << QString::number(vTop.y()) << " "
                << QString::number(vTop.z()) << "\n"
                << "    endloop\n  endfacet\n";

            normal = QVector3D::crossProduct(vNext-vCurrent, vCurrent-vBottom);
            out << "  facet normal "
                << QString::number(normal.x()) << " "
                << QString::number(normal.y()) << " "
                << QString::number(normal.z()) << "\n    outer loop\n"
                << "      vertex "
                << QString::number(vCurrent.x()) << " "
                << QString::number(vCurrent.y()) << " "
                << QString::number(vCurrent.z()) << "\n"
                << "      vertex "
                << QString::number(vBottom.x()) << " "
                << QString::number(vBottom.y()) << " "
                << QString::number(vBottom.z()) << "\n"
                << "      vertex "
                << QString::number(vNext.x()) << " "
                << QString::number(vNext.y()) << " "
                << QString::number(vNext.z()) << "\n"
                << "    endloop\n  endfacet\n";
        }
    }

    out << "endsolid Triangulation";

    file->close();
    delete file;
    if (m_consoleOutput)
        std::cout << "done." << std::endl;
}

void StlWriter::setAskOverwrite(bool newValue)
{
    m_askOverwrite = newValue;
}

bool StlWriter::getAskOverwrite()
{
    return m_askOverwrite;
}

void StlWriter::setConsoleOutput(bool newValue)
{
    m_consoleOutput = newValue;
}

bool StlWriter::getConsoleOutput()
{
    return m_consoleOutput;
}

void StlWriter::setDefaultFolder(QString defaultFolder)
{
    m_defaultFolder = defaultFolder;
}

QString StlWriter::getDefaultFolder()
{
    return m_defaultFolder;
}

void StlWriter::writeLineSegment(const QVector3D &start, const QVector3D &end, QTextStream *out, const double width)
{
    if (width == 0) {
        QVector3D v = (start+end)/2+QVector3D(0.05,0.05,0.05);

        QVector3D normal = QVector3D::crossProduct(end-start,v-start);
        normal.normalize();

        (*out) << "  facet normal "
               << QString::number(normal.x()) << " "
               << QString::number(normal.y()) << " "
               << QString::number(normal.z()) << "\n    outer loop\n"
               << "      vertex "
               << QString::number(start.x()) << " "
               << QString::number(start.y()) << " "
               << QString::number(start.z()) << "\n"
               << "      vertex "
               << QString::number(end.x()) << " "
               << QString::number(end.y()) << " "
               << QString::number(end.z()) << "\n"
               << "      vertex "
               << QString::number(v.x()) << " "
               << QString::number(v.y()) << " "
               << QString::number(v.z()) << "\n"
               << "    endloop\n  endfacet\n";
    } else {
        QVector3D normal1 = QVector3D::crossProduct(end-start, QVector3D(1,0,0));
        normal1.normalize();
        if (normal1.length() <= 0.1) {
            normal1 = QVector3D::crossProduct(end-start, QVector3D(0,1,0));
            normal1.normalize();
        }
        QVector3D normal2 = QVector3D::crossProduct(end-start, normal1);
        normal2.normalize();

        normal1 *= 0.5*width;
        normal2 *= 0.5*width;

        QVector3D vStart[4] = {start + normal1, start + normal2, start - normal1, start - normal2};
        QVector3D vEnd[4] = {end + normal1, end + normal2, end - normal1, end - normal2};
        QVector3D normal[4] = {(normal1+normal2).normalized(), (normal2-normal1).normalized(), (-normal1-normal2).normalized(), (normal1-normal2).normalized()};

        for (int i = 0; i < 4; i++) {
            int iNext = (i+1)%4;

            (*out) << "  facet normal "
                   << QString::number(normal[i].x()) << " "
                   << QString::number(normal[i].y()) << " "
                   << QString::number(normal[i].z()) << "\n    outer loop\n"
                   << "      vertex "
                   << QString::number(vStart[i].x()) << " "
                   << QString::number(vStart[i].y()) << " "
                   << QString::number(vStart[i].z()) << "\n"
                   << "      vertex "
                   << QString::number(vEnd[i].x()) << " "
                   << QString::number(vEnd[i].y()) << " "
                   << QString::number(vEnd[i].z()) << "\n"
                   << "      vertex "
                   << QString::number(vStart[iNext].x()) << " "
                   << QString::number(vStart[iNext].y()) << " "
                   << QString::number(vStart[iNext].z()) << "\n"
                   << "    endloop\n  endfacet\n";

            (*out) << "  facet normal "
                   << QString::number(normal[i].x()) << " "
                   << QString::number(normal[i].y()) << " "
                   << QString::number(normal[i].z()) << "\n    outer loop\n"
                   << "      vertex "
                   << QString::number(vEnd[i].x()) << " "
                   << QString::number(vEnd[i].y()) << " "
                   << QString::number(vEnd[i].z()) << "\n"
                   << "      vertex "
                   << QString::number(vEnd[iNext].x()) << " "
                   << QString::number(vEnd[iNext].y()) << " "
                   << QString::number(vEnd[iNext].z()) << "\n"
                   << "      vertex "
                   << QString::number(vStart[iNext].x()) << " "
                   << QString::number(vStart[iNext].y()) << " "
                   << QString::number(vStart[iNext].z()) << "\n"
                   << "    endloop\n  endfacet\n";
        }
    }
}

QFile *StlWriter::processFilename(QString filename)
{
    QString filepath = QString(m_defaultFolder).append(filename);
    QFile *file = new QFile(filepath);
    if (file->exists() && m_askOverwrite) {
        if (QMessageBox::question(0, QString("File already exists"),
                                  QString("The file %1 will be overwriten. Continue?").arg(filepath),
                                  QMessageBox::Yes | QMessageBox::No)
                == QMessageBox::No) {
            return 0;
        }
    }

    if (m_consoleOutput)
        std::cout << "Writing file:" << filepath.toStdString() << "...";

    return file;
}

bool StlWriter::removeRecursively(const QString &dirName)
{
    bool result = true;
    QDir dir(dirName);

    if (dir.exists(dirName)) {
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
            if (info.isDir()) {
                result = removeRecursively(info.absoluteFilePath());
            }
            else {
                result = QFile::remove(info.absoluteFilePath());
            }

            if (!result) {
                return result;
            }
        }
        result = dir.rmdir(dirName);
    }
    return result;
}

bool StlWriter::copyRecursively(const QString &srcFilePath, const QString &tgtFilePath)
{
    QFileInfo srcFileInfo(srcFilePath);
    if (srcFileInfo.isDir()) {
        QDir targetDir(tgtFilePath);
        targetDir.cdUp();
        if (!targetDir.mkpath(QFileInfo(tgtFilePath).fileName()))
            return false;
        QDir sourceDir(srcFilePath);
        QStringList fileNames = sourceDir.entryList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot | QDir::Hidden | QDir::System);
        foreach (const QString &fileName, fileNames) {
            const QString newSrcFilePath
                    = srcFilePath + QLatin1Char('/') + fileName;
            const QString newTgtFilePath
                    = tgtFilePath + QLatin1Char('/') + fileName;
            if (!copyRecursively(newSrcFilePath, newTgtFilePath))
                return false;
        }
    } else {
        if (!QFile::copy(srcFilePath, tgtFilePath))
            return false;
    }
    return true;
}

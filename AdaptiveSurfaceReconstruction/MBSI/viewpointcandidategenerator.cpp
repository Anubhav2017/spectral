#include "viewpointcandidategenerator.h"

#include "splineevaluator.h"

#include <cstdlib>
#include <ctime>

#include <QDebug>

ViewpointCandidateGenerator::ViewpointCandidateGenerator(BSplinePatchNetwork *bSplines):
    m_bSplinePatchNetwork(bSplines), m_cellMesh(bSplines->getCellMesh()), m_mesh(m_cellMesh->getOriginalMesh())
{
    srand(time(NULL));
}

void ViewpointCandidateGenerator::generateViewPoints_useVertices(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, const double cameraDistance)
{
    const int numberOfVertices = m_mesh->getNumberOfVertices();
    camPositions.resize(numberOfVertices);
    pivotPoints.resize(numberOfVertices);
    for (int iVert = 0; iVert < numberOfVertices; iVert++) {
        QVector3D vPos = m_mesh->getVertex(iVert)->getPosition();
        pivotPoints[iVert] = vPos;

        QVector3D normal = m_mesh->getVertexNormal(iVert).normalized();
        camPositions[iVert] = vPos + cameraDistance * normal;
    }
}

void ViewpointCandidateGenerator::generateViewPoints_random(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, const double cameraDistance, const int numberOfViewPoints, const bool randomPivotPoints)
{
    double minX, maxX, minY, maxY, minZ, maxZ;
    m_mesh->calculateMinMaxVertexPositions(minX, maxX, minY, maxY, minZ, maxZ);
    QVector3D center((float) (minX + maxX)/2, (float) (minY + maxY)/2, (float) (minZ + maxZ)/2);

    camPositions.resize(numberOfViewPoints);
    pivotPoints.resize(numberOfViewPoints);

    const int numberOfVertices = m_mesh->getNumberOfVertices();
    for (int i = 0; i < numberOfViewPoints; i++) {
        const float randX = -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2.0));
        const float randY = -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2.0));
        const float randZ = -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2.0));
        QVector3D randomDir(randX, randY, randZ);
        randomDir.normalize();

        if (randomPivotPoints)
            pivotPoints[i] = m_mesh->getVertex(rand() % numberOfVertices)->getPosition();
        else
            pivotPoints[i] = center;

        camPositions[i] = pivotPoints[i] + randomDir * cameraDistance;
    }
}

void ViewpointCandidateGenerator::generateViewPoints_subdivision(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, QVector<int> &depths, const double cameraDistance, const QVector<SubSurfaceTree *> subdivisions)
{
    foreach (SubSurfaceTree *subsurf, subdivisions) {
        BSplineSurface *surf = subsurf->getBaseSurface();
        QVector<SubSurfaceTree *> leaves = subsurf->getAllLeaves();
        foreach (SubSurfaceTree *leaf, leaves) {
            const double uCenter = leaf->getUCenter();
            const double vCenter = leaf->getVCenter();

            const QVector3D pivotPoint = surf->evaluate(uCenter, vCenter);
            const QVector3D normal = QVector3D::crossProduct(surf->evaluateTu(uCenter, vCenter), surf->evaluateTv(uCenter, vCenter)).normalized();
            const QVector3D camPosition = pivotPoint + cameraDistance * normal;

            //Add some random noise to position and/or radius
//            double uDiff = leaves[i]->getU1() - leaves[i]->getU0();
//            double vDiff = leaves[i]->getV1() - leaves[i]->getV0();
//            const double randomU = leaves[i]->getU0() + (((double) rand() / (RAND_MAX)) * uDiff);
//            const double randomV = leaves[i]->getV0() + (((double) rand() / (RAND_MAX)) * vDiff);
//            const double randomDistance = (0.95 + ((double) rand() / (RAND_MAX)) * 0.1) * cameraDistance;
//            const QVector3D pivotPoint = surf->evaluate(randomU, randomV);
//            const QVector3D normal = QVector3D::crossProduct(surf->evaluateTu(uCenter, vCenter), surf->evaluateTv(uCenter, vCenter)).normalized();
//            const QVector3D camPosition = pivotPoint + randomDistance * normal;

            camPositions << camPosition;
            pivotPoints << pivotPoint;
            depths << leaf->getDepth();
        }
    }
}

#ifndef VIEWPOINTCANDIDATEGENERATOR_H
#define VIEWPOINTCANDIDATEGENERATOR_H

#include "BSpline/bsplinepatchnetwork.h"
#include "subsurfacetree.h"

#include <QVector>
#include <QVector3D>
#include <QHash>
#include <QPair>

class ViewpointCandidateGenerator
{
public:
    ViewpointCandidateGenerator(BSplinePatchNetwork *bSplines);

    void generateViewPoints_useVertices(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, const double cameraDistance);
    void generateViewPoints_random(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, const double cameraDistance, const int numberOfViewPoints, const bool randomPivotPoints);

    void generateViewPoints_subdivision(QVector<QVector3D> &camPositions, QVector<QVector3D> &pivotPoints, QVector<int> &depths, const double cameraDistance, const QVector<SubSurfaceTree *> subdivisions);
private:
    BSplinePatchNetwork *m_bSplinePatchNetwork;
    CellMesh *m_cellMesh;
    Mesh *m_mesh;
};

#endif // VIEWPOINTCANDIDATEGENERATOR_H

#ifndef SURFACEVORONOIGENERATOR_H
#define SURFACEVORONOIGENERATOR_H

#include "CellMesh/cellmesh.h"
#include "Mesh/mesh.h"
#include "priorityqueue.h"

#include <QVector>
#include <QVector3D>

class SurfaceVoronoiGenerator
{
public:
    SurfaceVoronoiGenerator(Mesh *m, const bool initializeWithFirstFace = true);
    void initializeEmptyDijkstraVoronoi();
    CellMesh *generateCellMesh(bool atLeast3PlvPerEdge);

    void computeDiskLikeVoronoiTessellation(const int minNumberOfCells, bool instantTopologyCheck);

    bool isSource(Face *f);
    void addSource(Face *source);
    void addSourceWithLargestDistanceInCell(const int cellId);
    void addSourceFromLargestCellWithHigestDistance();
private:
    int calculateEulerNumberOfCell(const int cellId);
    int findLargestCell();
    int findHighestDijkstraDistIndexInCell(const int cellId);

    Mesh *m_mesh;
    const int m_numberOfMeshFaces;

    //Variables used by the algorithms
    QVector<Face *> m_sources;
    QVector<double> m_dijkstraDist;
    QVector<Face *> m_dijkstraPred;
    QVector<int> m_voronoiId;
    QVector<QVector<Face *> > m_voronoiCells;
    PriorityQueue<Face *> m_priorityQueue;
    QVector<bool> m_inQueue;

};

#endif // SURFACEVORONOIGENERATOR_H

#ifndef QUADMESHGENERATOR_H
#define QUADMESHGENERATOR_H

#include "Mesh/mesh.h"
#include "CellMesh/cellmesh.h"
#include "surfacevoronoigenerator.h"

class QuadmeshGenerator
{
public:
    QuadmeshGenerator(Mesh *m);

    CellMesh *computeQuadmeshByDelaunayMatching(bool instantTopologyCheckForVoronoi, int minNumberOfFinalCells, int minNumberOfVoronoiCells = 4);
    CellMesh *computeQuadmeshBySubdividingVoronoi(bool instantTopologyCheckForVoronoi, bool skipRefinement, int minNumberOfFinalCells, int minNumberOfVoronoiCells = 4);
    CellMesh *subdivideMeshByCuttingPlanes(QVector<QVector3D> &normals, QVector<QVector3D> &basePoints);

    void setParameterizationParameters(Parameterization::BorderType bt, Parameterization::WeightType wt, QPair<double, double> paramInterval, int storageReservation);
private:
    CellVertex *findCellVertexWithValencyUnequal3(CellMesh *cm);
    CellEdge *findDoubleCellEdgeBetweenCellFaces(CellMesh *cm);
    Face *findFreeVoronoiSourceAroundCellVertex(CellMesh *cm, CellVertex *cv, SurfaceVoronoiGenerator *voronoiGenerator);
    Face *findFreeVoronoiSourceAlongCellEdge(CellEdge *ce, SurfaceVoronoiGenerator *voronoiGenerator);
    void applyParametersToCellMesh(CellMesh *cm);

    Mesh *m_mesh;

    bool m_parameterizationParametersSet;
    Parameterization::BorderType m_paramBorderType;
    Parameterization::WeightType m_paramWeightType;
    QPair<double, double> m_paramParameterInterval;
    int m_paramStorageReservation;
};

#endif // QUADMESHGENERATOR_H

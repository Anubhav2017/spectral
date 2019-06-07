#ifndef CELLMESH_H
#define CELLMESH_H

#include "cellvertex.h"
#include "celledge.h"
#include "cellface.h"
#include "Mesh/mesh.h"
#include "parameterization.h"

#include <QVector>
#include <QList>
#include <QHash>

class CellMesh
{
public:
    CellMesh(Mesh *originalMesh, QPair<double, double> parameterInterval);
    CellMesh(Mesh *originalMesh, QVector<int> &faceCellMap, int numberOfCells, bool atLeast3PlvPerEdge = true, QPair<double, double> parameterInterval = QPair<double, double>(0,1));
    ~CellMesh();

    void straightenAllCellEdges();
    void optimizeAllVertices();
    void optimizeAllVerticesWithTemporarySubdivision(); //probably faster if almost all celledges are unsafe to be straightened directly
    CellEdge *straightenCellEdgeWithSafetyCheck(CellEdge *ce);
    CellEdge *straightenCellEdgeUnsafe(CellEdge *ce);
    CellEdge *straightenCellEdgeSafe(CellEdge *ce);
    CellVertex* optimizeVertexWithSafetyCheck(CellVertex *cv);
    CellVertex* optimizeVertexUnsafe(CellVertex *cv);
    CellVertex* optimizeVertexSafe(CellVertex *cv);
    void doSubdivisionForFullMesh();
    void doSubdivisionForOneCellFace(CellFace* cf);
    void transformToDualMesh();
    void computeFullParameterization(Parameterization::BorderType bt);

    //void cutCellBy3DPlane(const int id, const QVector3D point, const QVector3D normal);

    //Restrictions:
    //-The border of the input cell must not touch itself (i.e. each cell edge and vertex are contained only once in the cellface)
    //-The subdivided cells must have at least three cell vertices after the subdivision
    CellEdge *subdivideCellByCut(int id, PolylineVertex *vSeparationStart, PolylineVertex *vSeparationEnd);

    //Restriction:
    //-At least three separation vertices are required
    CellVertex *subdivideCellCentral(int id, QVector<PolylineVertex *> subdivisionVertices);
    CellFace *mergeCellsAroundEdge(CellEdge *separatingEdge);
    CellFace *mergeCellsAroundVertex(CellVertex *centerVertex);
    CellEdge *mergeEdgesAroundUnecessaryVertex(CellVertex *centerVertex);

    bool checkForConsistentOrientation();
    void checkAndEnforceCellOrientation();
    void enforceConsistentOrientation();
    void invertTotalOrientation();
    QVector<int> getValencyHistogramCellVertices();
    QVector<int> getValencyHistogramCellFaces();

    int getNumberOfVertices();
    CellVertex *getVertex(int id);

    int getNumberOfEdges();
    CellEdge *getEdge(int id);

    int getNumberOfFaces();
    CellFace *getFace(int id);

    int getNeighborId(CellFace *face, int number);
    CellFace *getNeighbor(CellFace *face, int number);

    void setParameterInterval(QPair<double, double> parameterInterval);
    QPair<double, double> getParameterInterval();

    Parameterization *getCellFaceParameterization(int i);
    Parameterization *getQuadCellParameterization(int i);
    QVector<double> getBoundaryParameterizationFromIncidentQuadParam(int ceId);
    void updateBoundaryParameterization(int ceId, const QVector<double> &newParam, const bool onlyInteriorPoints);
    void setParameterizationWeightType(Parameterization::WeightType wt);
    void setParameterizationBorderType(Parameterization::BorderType bt);
    int getParameterizationStorageReservation();
    void setParameterizationStorageReservation(int expectedEdgesPerVertex);
    double limitQuadCellParameterizationsToParameterIntervall();

    QVector<double> calculateApproximateSurfaceAreas();
    void translate(const QVector3D &translationVector);
    void scale(double factorX, double factorY, double factorZ);

    Mesh *getOriginalMesh();

    QString getInfo();
    void resetInfo();
    void addInfo(QString additionalInfo);

    void cleanUpMeshIndices();
    void checkMesh();
    void saveToFile(QString filename, QString optionalHeader = QString());
    static CellMesh *loadFromFile(QString cellmeshFile, Mesh *origMesh);
    void saveParameterizations(QString genericFilename, QString exportFilename = QString(), QString optionalHeader = QString());
    void loadParameterizations(QString genericFilename);
private:
    CellMesh(); //Constructor is private so it can only be accessed through static cell mesh methods
    int takeFreeVertexId();
    int takeFreeEdgeId();
    int takeFreeFaceId();
    CellVertex *subdivideEdge(CellEdge *ce, PolylineVertex *splitVertex, QVector<PolylineVertex *> *updateVertices = 0);
    QVector3D calculateIntersectionWithLine(Vertex *v0, Vertex *v1, Parameterization *param, const QVector2D &projectedStart, const QVector2D &normal, QVector2D *intersectionParameter);
    QVector2D calculateIntersectionCoefficients(const QVector2D &u0, const QVector2D &u1, const QVector2D &v0, const QVector2D &v1);
    void restrictToInterval(QVector2D *lambda, const float minX, const float maxX, const float minY, const float maxY);
    int findSlice(double angle, QVector<double> &separationAngles);

    QVector<CellVertex *> m_cellVertices;
    QVector<CellEdge *> m_cellEdges;
    QVector<CellFace *> m_cellFaces;
    QVector<Parameterization *> m_cellFaceParameterizations;
    QPair<double, double> m_parameterInterval;

    QList<int> m_freeVertexIds;
    QList<int> m_freeEdgeIds;
    QList<int> m_freeFaceIds;

    Mesh *m_originalMesh;

    Parameterization::WeightType m_parameterizationWeightType;
    Parameterization::BorderType m_parameterizationBorderType;
    int m_parameterizationStorageReservation_ExpectedEdgesPerVertex;

    QString m_info;
};

class PolylineVertexComparator
{
public:
    PolylineVertexComparator(QHash<PolylineVertex *, QVector2D> *param, QVector2D origin) :
        m_parameterization(param), m_origin(origin)
    {}

    bool operator () (PolylineVertex *plv0, PolylineVertex *plv1) {
        return ((m_parameterization->value(plv0) - m_origin).length() < (m_parameterization->value(plv1) - m_origin).length());
    }
private:
    QHash<PolylineVertex *, QVector2D> *m_parameterization;
    QVector2D m_origin;
};

#endif // CELLMESH_H

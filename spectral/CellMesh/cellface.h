#ifndef CELLFACE_H
#define CELLFACE_H

#include "Mesh/face.h"
#include "cellvertex.h"
#include "polylinevertex.h"

#include <QVector>

class CellFace : public Face
{
public:
    CellFace(QVector<Edge *> edges, QVector<Edge *> meshEdges, QVector<Vertex *> meshVertices, int id);

    QVector<Edge *> *getMeshEdges();
    QVector<Vertex *> *getMeshVertices();
    bool edgeIsInverted(int id);
    PolylineVertex *getPolyVertexCorrespondingToCellVertex(CellVertex *cv);

    void updateId(int newId);
    void replaceSplitEdges(Edge *old, Edge *firstNewEdge, Edge *secondNewEdge);
    void mergeEdges(Edge *oldEdge0, Edge *oldEdge1, Edge *newEdge);
    void switchOrientation();
private:
    void detectInvertedEdges();
    QVector<Edge *> m_meshEdges;
    QVector<Vertex *> m_meshVertices;
    QVector<bool> m_edgeInverted;
};

#endif // CELLFACE_H

#ifndef CELLEDGE_H
#define CELLEDGE_H

#include "Mesh/edge.h"
#include "cellvertex.h"
#include "polylinevertex.h"

class CellEdge : public Edge
{
public:
    CellEdge(CellVertex *v0, CellVertex *v1, QVector<PolylineVertex *> polylineVertices, int id);
    ~CellEdge();

    QVector<PolylineVertex *> *getPolylineVertices();
    QVector<QVector3D> getPolylineVertexPositions();

    void updateId(int newId);
    void updateIncidentFaceId(int oldId, int newId);
    void setG1Edge(bool newValue);
    bool isG1Edge();
private:
    QVector<PolylineVertex *> m_polylineVertices;   //polyline vertices are stored to begin at m_vertices[0] and end at m_vertices[1]
    bool m_G1Edge;
};

#endif // CELLEDGE_H

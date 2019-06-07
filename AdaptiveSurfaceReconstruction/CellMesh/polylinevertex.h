#ifndef POLYLINEVERTEX_H
#define POLYLINEVERTEX_H

#include "Mesh/vertex.h"
#include "Mesh/edge.h"

#include <QVector>
#include <QVector3D>

class PolylineVertex : public Vertex
{
public:
    PolylineVertex(QVector3D position, Edge *crossingEdge);
    PolylineVertex(Vertex *meshVertex);
    PolylineVertex(QVector3D position, int incidentFaceId);
    PolylineVertex(PolylineVertex *plv);

    bool isMeshVertex();
    bool isCrossingVertex();
    bool isFreeVertex();
    Vertex *getMeshVertex();
    Edge *getCrossingEdge();
private:
    Vertex *m_meshVertex;
    Edge *m_crossingEdge;
};

#endif // POLYLINEVERTEX_H

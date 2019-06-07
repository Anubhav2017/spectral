#include "polylinevertex.h"

PolylineVertex::PolylineVertex(QVector3D position, Edge *crossingEdge):
    Vertex(position, -1), m_meshVertex(0), m_crossingEdge(crossingEdge)
{
    m_incidentFaceIds << crossingEdge->getIncidentFaceId(0) << crossingEdge->getIncidentFaceId(1);
    m_incidentEdgeIds << crossingEdge->getId();
}

PolylineVertex::PolylineVertex(Vertex *meshVertex):
    Vertex(meshVertex->getPosition(), -1), m_meshVertex(meshVertex), m_crossingEdge(0)
{
    m_incidentFaceIds = *meshVertex->getIncidentFaceIds();
    m_incidentEdgeIds = *meshVertex->getIncidentEdgeIds();
}

PolylineVertex::PolylineVertex(QVector3D position, int incidentFaceId):
    Vertex(position, -1), m_meshVertex(0), m_crossingEdge(0)
{
    m_incidentFaceIds << incidentFaceId;
}

PolylineVertex::PolylineVertex(PolylineVertex *plv):
    Vertex(plv->getPosition(), plv->getId()), m_meshVertex(plv->m_meshVertex), m_crossingEdge(plv->m_crossingEdge)
{
    m_incidentFaceIds = *plv->getIncidentFaceIds();
    m_incidentEdgeIds = *plv->getIncidentEdgeIds();
}

bool PolylineVertex::isMeshVertex()
{
    return (m_meshVertex != 0);
}

bool PolylineVertex::isCrossingVertex()
{
    return (m_crossingEdge != 0);
}

bool PolylineVertex::isFreeVertex()
{
    return (m_meshVertex == 0 && m_crossingEdge == 0);
}

Vertex *PolylineVertex::getMeshVertex()
{
    return m_meshVertex;
}

Edge *PolylineVertex::getCrossingEdge()
{
    return m_crossingEdge;
}

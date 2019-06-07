#include "edge.h"

Edge::Edge(Vertex *v0, Vertex *v1, int id):
    m_id(id)
{
    m_vertices[0] = v0;
    m_vertices[1] = v1;

    m_incidentFaceIds[0] = -1;
    m_incidentFaceIds[1] = -1;
}

int Edge::getId()
{
    return m_id;
}

Vertex *Edge::getVertex(int id)
{
    return m_vertices[id];
}

Vertex *Edge::getConnectedVertex(Vertex *v)
{
    if (m_vertices[0] == v)
        return m_vertices[1];
    else if (m_vertices[1] == v)
        return m_vertices[0];

    return 0;
}

int Edge::getConnectedFaceId(int fId)
{
    if (m_incidentFaceIds[0] == fId)
        return m_incidentFaceIds[1];
    else
        return m_incidentFaceIds[0];
}

bool Edge::isIncident(Vertex *v)
{
    return (m_vertices[0] == v || m_vertices[1] == v);
}

bool Edge::addIncidentFaceId(int id)
{
    if (m_incidentFaceIds[0] == -1)
        m_incidentFaceIds[0] = id;
    else if (m_incidentFaceIds[1] == -1)
        m_incidentFaceIds[1] = id;
    else
        return false;
    return true;
}

int Edge::getIncidentFaceId(int number)
{
    return m_incidentFaceIds[number];
}

Vertex *Edge::getSharedVertex(Edge *neighbor)
{
    if (neighbor->getVertex(0) == m_vertices[0] || neighbor->getVertex(1) == m_vertices[0])
        return m_vertices[0];
    else if (neighbor->getVertex(0) == m_vertices[1] || neighbor->getVertex(1) == m_vertices[1])
        return m_vertices[1];

    return 0;
}

Vertex *Edge::getUnsharedVertex(Edge *neighbor)
{
    if (neighbor->getVertex(0) == m_vertices[0] || neighbor->getVertex(1) == m_vertices[0])
        return m_vertices[1];
    else if (neighbor->getVertex(0) == m_vertices[1] || neighbor->getVertex(1) == m_vertices[1])
        return m_vertices[0];

    return 0;
}

double Edge::getDistanceBetweenVertices()
{
    return (m_vertices[0]->getPosition() - m_vertices[1]->getPosition()).length();
}


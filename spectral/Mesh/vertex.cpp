#include "vertex.h"

Vertex::Vertex(QVector3D position, int id):
    m_id(id), m_pos(position)
{
}

int Vertex::getId()
{
    return m_id;
}

float Vertex::getPosX()
{
    return m_pos.x();
}

float Vertex::getPosY()
{
    return m_pos.y();
}

float Vertex::getPosZ()
{
    return m_pos.z();
}

QVector3D Vertex::getPosition()
{
    return m_pos;
}

void Vertex::updatePosition(QVector3D newPosition)
{
    m_pos = newPosition;
}

int Vertex::getConnectingEdgeId(Vertex *otherVertex)
{
    foreach (int eId, *otherVertex->getIncidentEdgeIds())
        if (m_incidentEdgeIds.contains(eId))
            return eId;

    return -1;
}

int Vertex::getSharedFaceId(Vertex *otherVertex)
{
    foreach (int fId, *otherVertex->getIncidentFaceIds())
        if (m_incidentFaceIds.contains(fId))
            return fId;

    return -1;
}

void Vertex::addIncidentEdgeId(int id)
{
    m_incidentEdgeIds.append(id);
}

QVector<int> *Vertex::getIncidentEdgeIds()
{
    return &m_incidentEdgeIds;
}

void Vertex::addIncidentFaceId(int id)
{
    m_incidentFaceIds.append(id);
}

QVector<int> *Vertex::getIncidentFaceIds()
{
    return &m_incidentFaceIds;
}


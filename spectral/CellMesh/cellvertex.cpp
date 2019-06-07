#include "cellvertex.h"

CellVertex::CellVertex(QVector3D position, int id) :
    Vertex(position, id)
{
}

void CellVertex::updateId(int newId)
{
    m_id = newId;
}

void CellVertex::updateIncidentEdgeId(int oldId, int newId)
{
    const int index = m_incidentEdgeIds.indexOf(oldId);
    if (index != -1) {
        if (newId != -1) {
            m_incidentEdgeIds[index] = newId;
        } else {
            m_incidentEdgeIds.remove(index);
        }
    }
}

void CellVertex::updateIncidentFaceId(int oldId, int newId)
{
    const int index = m_incidentFaceIds.indexOf(oldId);
    if (index != -1) {
        if (newId != -1) {
            m_incidentFaceIds[index] = newId;
        } else {
            m_incidentFaceIds.remove(index);
        }
    }
}

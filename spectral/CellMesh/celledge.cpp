#include "celledge.h"

CellEdge::CellEdge(CellVertex *v0, CellVertex *v1, QVector<PolylineVertex *> polylineVertices, int id) :
    Edge(v0, v1, id), m_polylineVertices(polylineVertices), m_G1Edge(true)
{
    PolylineVertex *startVertex = polylineVertices.first();
    PolylineVertex *endVertex = polylineVertices.last();

//    if (polylineVertices.size() == 2) {
//        qDebug() << "inserted";
//        polylineVertices.insert(1, new PolylineVertex((startVertex->getPosition() + endVertex->getPosition())/2, -1));
//    }

    if (startVertex->getPosition() == v1->getPosition() && endVertex->getPosition() == v0->getPosition()) {
        m_vertices[0] = v1;
        m_vertices[1] = v0;
    }

//    if (startVertex->getPosition() != v0->getPosition() || endVertex->getPosition() != v1->getPosition()) {
//        qDebug() << "ERROR IN EDGE CREATION";
//        qDebug() << startVertex->getPosition() << endVertex->getPosition() << v0->getPosition() << v1->getPosition();
//    }
//    qDebug() << "Debug in CellEdge constructor";
}

CellEdge::~CellEdge()
{
    foreach (PolylineVertex *plv, m_polylineVertices)
        delete plv;
}

QVector<PolylineVertex *> *CellEdge::getPolylineVertices()
{
    return &m_polylineVertices;
}

QVector<QVector3D> CellEdge::getPolylineVertexPositions()
{
    const int numPLV = m_polylineVertices.size();
    QVector<QVector3D> positions(numPLV);
    for (int i = 0; i < numPLV; i++)
        positions[i] = m_polylineVertices[i]->getPosition();

    return positions;
}

void CellEdge::updateId(int newId)
{
    m_id = newId;
}

void CellEdge::updateIncidentFaceId(int oldId, int newId)
{
    if (m_incidentFaceIds[0] == oldId)
        m_incidentFaceIds[0] = newId;
    else if (m_incidentFaceIds[1] == oldId)
        m_incidentFaceIds[1] = newId;
}

void CellEdge::setG1Edge(bool newValue)
{
    m_G1Edge = newValue;
}

bool CellEdge::isG1Edge()
{
    return m_G1Edge;
}

#include "cellface.h"
#include "celledge.h"

#include <QDebug>

CellFace::CellFace(QVector<Edge *> edges, QVector<Edge *> meshEdges, QVector<Vertex *> meshVertices, int id) :
    Face(edges, QVector3D(0,0,0), id), m_meshEdges(meshEdges), m_meshVertices(meshVertices), m_edgeInverted(QVector<bool>(edges.size(), false))
{
    this->extractVerticesFromEdges();
    this->detectInvertedEdges();
}

QVector<Edge *> *CellFace::getMeshEdges()
{
    return &m_meshEdges;
}

QVector<Vertex *> *CellFace::getMeshVertices()
{
    return &m_meshVertices;
}

bool CellFace::edgeIsInverted(int id)
{
    return m_edgeInverted[id];
}

PolylineVertex *CellFace::getPolyVertexCorrespondingToCellVertex(CellVertex *cv)
{
    //Edge i connects vertices i and i+1 -> if it is not inverted, return the last polylinevertex
    const int index = m_vertices.indexOf(cv);
    CellEdge *ce = (CellEdge *) m_edges[index];
    if (edgeIsInverted(index)) {
        return ce->getPolylineVertices()->last();
    } else {
        return ce->getPolylineVertices()->first();
    }
}

void CellFace::updateId(int newId)
{
    m_id = newId;
}

void CellFace::replaceSplitEdges(Edge *old, Edge *firstNewEdge, Edge *secondNewEdge)
{
    const int index = m_edges.indexOf(old);
    m_edges[index] = secondNewEdge;     //First, insert the second edge (will be moved back, when the first edge is inserted)
    m_edges.insert(index, firstNewEdge);

    if (!m_edges[(index+2) % m_edges.size()]->getSharedVertex(secondNewEdge)) {
        m_edges[index] = secondNewEdge;
        m_edges[index+1] = firstNewEdge;
    }

    //Not the most efficient solution, but probably the most bug-free option
    m_vertices.insert(index+1, firstNewEdge->getSharedVertex(secondNewEdge));
    this->detectInvertedEdges();
}

void CellFace::mergeEdges(Edge *oldEdge0, Edge *oldEdge1, Edge *newEdge)
{
    m_edges[m_edges.indexOf(oldEdge0)] = newEdge;
    m_edges.remove(m_edges.indexOf(oldEdge1));

    this->extractVerticesFromEdges();
    this->detectInvertedEdges();
}

void CellFace::switchOrientation()
{
    const int numEdges = m_edges.size();
    QVector<Edge *> invertedEdges(numEdges);
    for (int i = 0; i < numEdges; i++)
        invertedEdges[i] = m_edges[numEdges - 1 - i];

    m_edges = invertedEdges;

    this->extractVerticesFromEdges();
    this->detectInvertedEdges();
}

void CellFace::detectInvertedEdges()
{
    m_edgeInverted = QVector<bool>(m_edges.size());
    const int numberOfVertices = m_vertices.size();
    for (int i = 0; i < m_edges.size(); i++) {
        CellEdge *ce = (CellEdge *) m_edges[i];
        if (ce->getVertex(0)->getPosition() == m_vertices[i]->getPosition() &&
                ce->getVertex(1)->getPosition() == m_vertices[(i+1)%numberOfVertices]->getPosition())
            m_edgeInverted[i] = false;
        else if (ce->getVertex(1)->getPosition() == m_vertices[i]->getPosition() &&
                 ce->getVertex(0)->getPosition() == m_vertices[(i+1)%numberOfVertices]->getPosition())
            m_edgeInverted[i] = true;
//        else
//            qDebug() << "Something went wrong with edge inversion detection (in CellFace)" << m_id << ce->getVertex(0)->getPosition() << ce->getVertex(1)->getPosition() << m_vertices[i]->getPosition() << m_vertices[(i+1)%numberOfVertices]->getPosition();
    }

}

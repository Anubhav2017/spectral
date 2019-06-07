#include "face.h"

Face::Face(QVector<Edge *> edges, QVector3D normal, int id):
    m_id(id), m_normal(normal), m_edges(edges)
{
    this->extractVerticesFromEdges();
}

int Face::getId()
{
    return m_id;
}

int Face::getNumberOfEdges()
{
    return m_edges.size();
}


Edge *Face::getEdge(int id)
{
    return m_edges[id];
}

int Face::getInternalEdgeId(Edge *e)
{
    return m_edges.indexOf(e);
}

int Face::getNumberOfVertices()
{
    return m_vertices.size();
}

Vertex *Face::getVertex(int id)
{
    return m_vertices[id];
}

int Face::getInternalVertexId(Vertex *v)
{
    return m_vertices.indexOf(v);
}

int Face::getLocalEdgeId(Edge *e)
{
    return m_edges.indexOf(e);
}

float Face::getNormalX()
{
    return m_normal.x();
}

float Face::getNormalY()
{
    return m_normal.y();
}

float Face::getNormalZ()
{
    return m_normal.z();
}

QVector3D Face::getNormal()
{
    return m_normal;
}

void Face::recomputeNormal()
{
    QVector3D newNormal = QVector3D::crossProduct(m_vertices[0]->getPosition() - m_vertices[1]->getPosition(), m_vertices[2]->getPosition() - m_vertices[1]->getPosition());
    if (QVector3D::dotProduct(newNormal, m_normal) < 0)
        newNormal *= -1;
    m_normal = newNormal.normalized();
}

QVector3D Face::getCenter()
{
    QVector3D center;
    foreach (Vertex *v, m_vertices)
        center += v->getPosition();

    center /= m_vertices.size();
    return center;
}

double Face::calculateAreaByCorners()
{
    const int numCorners = m_vertices.size();
    QVector3D corner0 = m_vertices[0]->getPosition();
    double result = 0;
    for (int i = 2; i < numCorners; i++)
        result += QVector3D::crossProduct(corner0 - m_vertices[i - 1]->getPosition(), corner0 - m_vertices[i]->getPosition()).length()/2;

    return result;
}

Edge *Face::getSharedEdge(Face *neighbor)
{
    foreach (Edge *e, m_edges)
        if (neighbor->getInternalEdgeId(e) != -1)
            return e;
    return 0;
}

void Face::extractVerticesFromEdges()
{
    const int numberOfEdges = m_edges.size();
    m_vertices.resize(numberOfEdges);

    if (numberOfEdges == 2) {
        m_vertices[0] = m_edges[0]->getVertex(0);
        m_vertices[1] = m_edges[0]->getVertex(1);
    } else if (numberOfEdges > 2){
        m_vertices[0] = m_edges.first()->getSharedVertex(m_edges.last());
        for (int i = 1; i < numberOfEdges; i++)
            m_vertices[i] = m_edges[i]->getSharedVertex(m_edges[i-1]);
    } else if (numberOfEdges == 1) {
        m_vertices[0] = m_edges[0]->getVertex(0);
        if (m_edges[0]->getVertex(1) != m_vertices[0])
            m_vertices.append(m_edges[0]->getVertex(1));
    }   //only case left is numberOfEdges == 0 -> nothing to do then
}


#ifndef FACE_H
#define FACE_H

#include "edge.h"
#include "vertex.h"

#include <QList>
#include <QVector3D>

class Face
{
public:
    Face(QVector<Edge *> edges, QVector3D normal, int id);

    int getId();

    int getNumberOfEdges();
    Edge *getEdge(int id);
    int getInternalEdgeId(Edge *e);

    int getNumberOfVertices();
    Vertex *getVertex(int id);
    int getInternalVertexId(Vertex *v);

    int getLocalEdgeId(Edge *e);

    float getNormalX();
    float getNormalY();
    float getNormalZ();
    QVector3D getNormal();
    void recomputeNormal();

    QVector3D getCenter();

    double calculateAreaByCorners();

    Edge *getSharedEdge(Face *neighbor);
protected:
    void extractVerticesFromEdges();
    int m_id;
    QVector3D m_normal;

    QVector<Edge *> m_edges;        //Edge i connects vertices i and i+1
    QVector<Vertex *> m_vertices;   //Vertex i is incident to edges i and i-1

};

#endif // FACE_H

#ifndef EDGE_H
#define EDGE_H

#include "vertex.h"

class Edge
{
public:
    Edge(Vertex *v0, Vertex *v1, int id);

    int getId();

    Vertex *getVertex(int id);
    Vertex *getConnectedVertex(Vertex *v);
    int getConnectedFaceId(int fId);
    bool isIncident(Vertex *v);

    bool addIncidentFaceId(int id);
    int getIncidentFaceId(int number);

    Vertex *getSharedVertex(Edge *other);
    Vertex *getUnsharedVertex(Edge *neighbor);

    double getDistanceBetweenVertices();
protected:
    int m_id;
    Vertex *m_vertices[2];
    int m_incidentFaceIds[2];
};

#endif // EDGE_H

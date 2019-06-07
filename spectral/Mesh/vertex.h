#ifndef VERTEX_H
#define VERTEX_H

#include <QVector>
#include <QVector3D>

class Vertex
{
public:
    Vertex(QVector3D position, int id);

    int getId();

    float getPosX();
    float getPosY();
    float getPosZ();
    QVector3D getPosition();
    void updatePosition(QVector3D newPosition);

    //Returns the id of the edge/face connecting this vertex to otherVertex
    //Returns -1 if no such edge/face exists
    //Returns the first found ID if several edges/faces exist
    int getConnectingEdgeId(Vertex *otherVertex);
    int getSharedFaceId(Vertex *otherVertex);

    void addIncidentEdgeId(int id);
    QVector<int> *getIncidentEdgeIds();
    void addIncidentFaceId(int id);
    QVector<int> *getIncidentFaceIds();
protected:
    int m_id;
    QVector3D m_pos;

    QVector<int> m_incidentEdgeIds;
    QVector<int> m_incidentFaceIds;
};

#endif // VERTEX_H

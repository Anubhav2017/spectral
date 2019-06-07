#ifndef CELLVERTEX_H
#define CELLVERTEX_H

#include "Mesh/vertex.h"

class CellVertex : public Vertex
{
public:
    CellVertex(QVector3D position, int id);

    void updateId(int newId);
    void updateIncidentEdgeId(int oldId, int newId);
    void updateIncidentFaceId(int oldId, int newId);
};

#endif // CELLVERTEX_H

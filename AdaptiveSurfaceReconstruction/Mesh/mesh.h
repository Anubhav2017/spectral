#ifndef MESH_H
#define MESH_H

#include "vertex.h"
#include "edge.h"
#include "face.h"

#include <QVector>
#include <QHash>
#include <QVector3D>
#include <QPair>
#include <QString>

class Mesh
{
public:
    Mesh();
    ~Mesh();

    Face *addFace(QVector<QVector3D> corners, QVector3D normal);

    int getNumberOfVertices();
    Vertex *getVertex(int id);

    int getNumberOfEdges();
    Edge *getEdge(int id);

    int getNumberOfFaces();
    Face *getFace(int id);

    Face *getNeighbor(Face *f, int id);
    QVector3D getVertexNormal(int vertexId);

    QString getInfo();
    void resetInfo();
    void addInfo(QString additionalInfo);

    int getMaxNumberOfEdgesOnVertex();
    void calculateMinMaxVertexPositions(double &minX, double &maxX, double &minY, double &maxY, double &minZ, double &maxZ);
    void calculateMinimumBoundingSphere(double &radius, QVector3D &center);
    void scaleMeshToFitIntoCube(bool keepProportions = true, double cubeEdgeLength = 1);
    void scaleMeshToFitIntoSphere(double radius = 1, QVector3D center = QVector3D(0, 0, 0));

    double calculateTotalSurfaceArea();
    void smoothSurface(const int maxIterations, const double maxDistFromOrigin, const double abortTotalDistance);

    void saveToFile(QString filename, QString optionalHeader = QString());
    static Mesh *loadFromFile(QString filename);
    static Mesh *loadFromBinaryStl(QString filename);

private:
    QHash<QVector3D, Vertex *> m_coordinateVertexMap;
    QHash<QPair<Vertex *, Vertex *>, Edge *> m_vertexPairEdgeMap;

    QVector<Vertex *> m_vertices;
    QVector<Edge *> m_edges;
    QVector<Face *> m_faces;

    QString m_info;
};

#endif // MESH_H

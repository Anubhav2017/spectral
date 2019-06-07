#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H
#include<Mesh/vertex.h>
#include<Mesh/edge.h>
#include<Mesh/mesh.h>
#include<Mesh/face.h>
class LaplaceSolver
{
    QVector<QVector<double>*> Laplacian;
    QVector<QVector<double>*> weights;
    int N;
    int numberOfEdges;
    QVector<Edge*> edges;
public:
    LaplaceSolver(Mesh* m);
    double LaplaceSolver::getWeights(Edge *e, Vertex *optionalOuterTriangleVertex);



};

#endif // LAPLACESOLVER_H

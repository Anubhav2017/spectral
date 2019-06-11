#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H
#include<Mesh/vertex.h>
#include<Mesh/edge.h>
#include<Mesh/mesh.h>
#include<Mesh/face.h>
#include<Eigen/Eigenvalues>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<iostream>

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
typedef Matrix<std::complex<double>,Dynamic,1> VectorXd;

class LaplaceSolver
{
    MatrixXd Laplacian;
    QVector<QVector<double> > weights;
    int N;
    int numberOfEdges;
    Mesh* m_originalMesh;
    QVector<Edge*> edges;
    QVector<Eigen::VectorXcd> functions;
    Eigen::VectorXcd eigenvalues;
public:
    LaplaceSolver(Mesh* m);
    double getWeights(Edge *e, Vertex *optionalOuterTriangleVertex);
    void decompose();
    Eigen::VectorXcd* getEigenValues();
    Eigen::VectorXcd* getEigenFunction(int id);
    std::complex<double> getEigenValue(int id);



};

#endif // LAPLACESOLVER_H

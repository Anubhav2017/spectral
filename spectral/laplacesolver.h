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
#include <QColor>
#include <QString>
#include <QFile>
#include <QTextStream>

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
    QVector<QColor> m_colorMapColors;
    QVector<double> m_colorMapPercentiles;

public:
    LaplaceSolver(Mesh* m);
    double getWeights(Edge *e, Vertex *optionalOuterTriangleVertex);
    void decompose();
    const Eigen::VectorXcd* getEigenValues() const;
    const Eigen::VectorXcd* getEigenFunction(int id) const;
    const std::complex<double> getEigenValue(int id) const;
    const QVector<QColor> generateColorMap(int id) const;
    const QVector<QColor> generateColorMap2(int id) const;

    void writeMeshWithVertexColors(const QVector<QColor> &colormap, const QString filename);



};

#endif // LAPLACESOLVER_H

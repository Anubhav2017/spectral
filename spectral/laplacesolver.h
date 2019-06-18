#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H
#include<Mesh/vertex.h>
#include<Mesh/edge.h>
#include<Mesh/mesh.h>
#include<Mesh/face.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/SparseCore>
#include<Eigen/Sparse>
#include<iostream>
#include <QColor>
#include <QString>
#include <QFile>
#include <QTextStream>
#include<QDebug>
#include<Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include<Spectra/MatOp/SparseGenMatProd.h>
//typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> VectorXd;
typedef Eigen::Triplet<double> T;

using namespace Spectra;

class LaplaceSolver
{
    QVector<T> tripletList;
    Eigen::SparseMatrix<double> Laplacian;
    QVector<QVector<double> > weights;
    int N;
    int numberOfEdges;
    Mesh* m_originalMesh;
    QVector<Edge*> edges;
    Eigen::MatrixXcd functions;
    Eigen::VectorXcd eigenvalues;
    QVector<QColor> m_colorMapColors;
    QVector<double> m_colorMapPercentiles;
    QVector<double> diagonalentries;

public:
    LaplaceSolver(Mesh* m);
    double getWeights(Edge *e, Vertex *optionalOuterTriangleVertex);
    void decompose();
    const Eigen::VectorXcd* getEigenValues() const;
    const Eigen::VectorXcd getEigenFunction(int id) const;
    const std::complex<double> getEigenValue(int id) const;
    const QVector<QColor> generateColorMap(int id) const;
    const QVector<QColor> generateColorMap2(int id) const;

    void writeMeshWithVertexColors(const QVector<QColor> &colormap, const QString filename);
    void writeVectors();
    const QVector<QColor> generateColorMap2(QVector<std::complex<double> >eigenvector) const;
    const void writeVector(Eigen::VectorXcd v, QString filename)const;
    const void writeMeshWithVector(const Eigen::VectorXcd v, const QString filename) const;



};

#endif // LAPLACESOLVER_H

#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include <QHash>
#include <QVector>
#include <QVector2D>
#include <QString>

#include "Mesh/mesh.h"
#include "Mesh/vertex.h"
#include "Mesh/edge.h"
#include "Mesh/face.h"
#include "CellMesh/cellface.h"
#include "CellMesh/polylinevertex.h"

#include "IterativeLinearSolvers"   //includes least squares solver which is only dev. version right now
#include <SparseLU>

class Parameterization
{
public:
    enum WeightType {WT_constant = 0, WT_length = 1, WT_cotangent = 2};
    enum BorderType {BT_equidistantCircle = 0, BT_lengthWeightedCircle = 1, BT_equidistantPolygon = 2, BT_equidistantUnitSquare = 3, BT_lengthWeightedUnitSquare = 4};

    //TODO add lengthWeightedUnitSquare

    bool contains(Vertex *v);
    void updateParameter(Vertex *v, double newU, double newV);
    QVector2D getParameter(Vertex *v);
    QVector<double> getParameterizationOfBoundaryCurve(int id);
    QVector2D getPolarCoordinates(Vertex *v);
    double getPolarAngle(Vertex *v);
    void saveAsImage(QString filename, int size = 1000, int symbolRadius1 = 2, int symbolRadius2 = 4);
    void addCopy(Vertex *original, Vertex *copy);
    Parameterization::BorderType getBorderType();

    void saveToFile(QString filename, QString optionalHeader = QString());
    static Parameterization *loadFromFile(QString filename, CellFace *cf, Mesh *originalMesh);

    QPair<double, double> getParameterInterval();
    QString getInfo();
    void addInfo(QString info);

    double limitToParameterInterval();
    double limitTo(const double lowerBoundU = 0, const double upperBoundU = 1, const double lowerBoundV = 0, const double upperBoundV = 1);

    static Parameterization *harmonicMap(CellFace *cf, Mesh *originalMesh, Parameterization::WeightType wt, Parameterization::BorderType bt, QPair<double, double> parameterInterval, const int expectedEdgesPerVertex = 12);
private:
    Parameterization();
    double getWeight(Edge *e, Parameterization::WeightType wt, Vertex *optionalOuterTriangleVertex = 0);
    static void calculateBorderProjection(QVector<PolylineVertex *> *borderVertices, QVector<int> *corners, Eigen::VectorXd &projectedX, Eigen::VectorXd &projectedY, Parameterization::BorderType bt, QPair<double, double> parameterInterval);
    void fillVertexLocalIdMap();

    QString m_info;
    QPair<double, double> m_parameterInterval;
    QHash<Vertex *, int> m_vertexLocalIdMap;
    Eigen::VectorXd m_parameterValuesU;
    Eigen::VectorXd m_parameterValuesV;
    CellFace *m_cellFace;
    Mesh *m_originalMesh;
    int m_numberOfVertices;
    int m_numberOfBorderVertices;
    int m_numberOfInnerVertices;
    Parameterization::BorderType m_borderType;
};

#endif // PARAMETERIZATION_H

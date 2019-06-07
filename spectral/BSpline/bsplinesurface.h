#ifndef BSPLINESURFACE_H
#define BSPLINESURFACE_H

#include <QVector>
#include <QVector3D>
#include <QPair>

#include "bsplineproperties.h"

class BSplineSurface
{
public:
    BSplineSurface(const BSplineProperties * const propertiesU, const BSplineProperties * const propertiesV, const QVector<QVector<QVector3D> > &controlPoints);

    QVector3D evaluate(double u, double v); //TODO implement evaluation with DeBoor Algorithm (should be way faster)
    QVector3D evaluateTu(double u, double v);
    QVector3D evaluateTv(double u, double v);
    QVector3D evaluateDerivative(double u, double v, int du, int dv);
    void evaluateFirstFundamentalForm(double u, double v, double &E, double &F, double &G);
    void evaluateSecondFundamentalForm(double u, double v, double &L, double &M, double &N);
    double evaluateGaussianCurvature(double u, double v);
    double evaluateMeanCurvature(double u, double v);
    void evaluateGaussianAndMeanCurvature(double u, double v, double &K, double &H);
    void evaluatePrincipalCurvatures(double u, double v, double &k1, double &k2);
    void translate(const QVector3D translationVector);
    void scale(double factorX, double factorY, double factorZ);

    BSplineProperties const * getPropertiesU() const;
    BSplineProperties const * getPropertiesV() const;
    QVector<QVector<QVector3D> > *getControlPointsP();
    QVector<QVector3D> getControlPointsAsVector();

    static int matToVecIndexLocal(const int &i, const int &j, const int &numberOfControlPointsPerRow);
    static int matToVecIndexLocal(const QPair<int, int> &ij, const int &numberOfControlPointsPerRow);
    static void vecToMatIndexLocal(const int &iVec, int &iOut, int &jOut, const int &numberOfControlPointsPerRow);
private:
    BSplineProperties const * const m_propertiesU;
    BSplineProperties const * const m_propertiesV;
    QVector<QVector<QVector3D> > m_controlPoints;
};

#endif // BSPLINESURFACE_H

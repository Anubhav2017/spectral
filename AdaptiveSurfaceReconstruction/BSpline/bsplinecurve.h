#ifndef BSPLINE_H
#define BSPLINE_H

#include <QVector3D>
#include <QVector>

#include "bsplineproperties.h"

class BSplineCurve
{
public:
    BSplineCurve(const BSplineProperties * const properties, const QVector<QVector3D> &controlPoints);

    void translate(const QVector3D &translationVector);
    void scale(const double factorX, const double factorY, const double factorZ);

    QVector3D evaluate(const double x);
    QVector3D evaluateTangent(const double x);
    QVector3D evaluateAnalytical(const double x);
    QVector3D evaluateDeBoor(const double x);

    double calculateLengthNumerically(const double x0, const double x1, const int numberOfSamplePoints = 1001);

    const BSplineProperties* getProperties() const;
    QVector<QVector3D> getControlPoints() const;
    const QVector<QVector3D>* getControlPointsP() const;
private:
    QVector3D d(const int i, const int k, const double x);

    const BSplineProperties * const m_properties;
    QVector<QVector3D> m_controlPoints;
};

#endif // BSPLINE_H

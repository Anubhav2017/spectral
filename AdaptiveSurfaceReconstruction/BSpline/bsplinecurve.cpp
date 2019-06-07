#include "bsplinecurve.h"

#include <QDebug>

BSplineCurve::BSplineCurve(const BSplineProperties * const properties, const QVector<QVector3D> &controlPoints) :
    m_properties(properties), m_controlPoints(controlPoints)
{

}

void BSplineCurve::translate(const QVector3D &translationVector)
{
    const int numCtrPoints = m_controlPoints.size();
    for (int i = 0; i < numCtrPoints; i++)
        m_controlPoints[i] += translationVector;
}

void BSplineCurve::scale(const double factorX, const double factorY, const double factorZ)
{
    const int numCtrPoints = m_controlPoints.size();
    for (int i = 0; i < numCtrPoints; i++) {
        QVector3D pos = m_controlPoints[i];
        pos.setX(pos.x() * factorX);
        pos.setY(pos.y() * factorY);
        pos.setZ(pos.z() * factorZ);
        m_controlPoints[i] = pos;
    }
}

QVector3D BSplineCurve::evaluate(const double x)
{
    return this->evaluateDeBoor(x);
}

QVector3D BSplineCurve::evaluateTangent(const double x)
{
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++) {
        result += m_controlPoints[i] * m_properties->NDerivative(x, i, 1);
    }
    return result;
}

QVector3D BSplineCurve::evaluateAnalytical(const double x)
{
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++) {
        result += m_controlPoints[i] * m_properties->N(x, i);
    }
    return result;
}

QVector3D BSplineCurve::evaluateDeBoor(const double x)
{
    //find l
    int l = 0;
    const QVector<double> *knots = m_properties->getKnotsP();
    while (l < m_controlPoints.size()-1 && !(x >= knots->at(l) && x < knots->at(l+1))) {
        l++;
    }

    return d(l, m_properties->getOrder(), x);
}

double BSplineCurve::calculateLengthNumerically(const double x0, const double x1, const int numberOfSamplePoints)
{
    QVector3D lastCoordinate = this->evaluate(x0);
    double result = 0;
    for (int i = 1; i < numberOfSamplePoints; i++) {
        const double factor = (double) i / (numberOfSamplePoints - 1);
        const QVector3D currentCoordinate = this->evaluate(factor * x1 + (1 - factor) * x0);
        result += (currentCoordinate - lastCoordinate).length();
        lastCoordinate = currentCoordinate;
    }

    return result;
}

const BSplineProperties *BSplineCurve::getProperties() const
{
    return m_properties;
}

QVector<QVector3D> BSplineCurve::getControlPoints() const
{
    return m_controlPoints;
}

const QVector<QVector3D> *BSplineCurve::getControlPointsP() const
{
    return &m_controlPoints;
}

QVector3D BSplineCurve::d(const int i, const int k, const double x)
{
    if (k == 0) {
        return m_controlPoints[i];
    } else {
        const QVector<double> *knots = m_properties->getKnotsP();
        double alpha = (x - knots->at(i)) / (knots->at(i + m_properties->getOrder() + 1 - k) - knots->at(i));
        return (d(i - 1, k - 1, x) * (1 - alpha)) + (d(i, k - 1, x) * alpha);
    }
}

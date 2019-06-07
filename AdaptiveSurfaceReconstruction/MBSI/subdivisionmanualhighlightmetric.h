#ifndef SUBDIVISIONMANUALHIGHLIGHTMETRIC_H
#define SUBDIVISIONMANUALHIGHLIGHTMETRIC_H

#include "subdivisionabstractmetric.h"

class SubdivisionManualHighlightMetric : public SubdivisionAbstractMetric
{
public:
    SubdivisionManualHighlightMetric(const int sampleRatePerKnotInterval);

    void addHighlightBoundingBox(const QVector3D &corner0, const QVector3D &corner1, const double highlightValue);
    void addHighlightBoundingSphere(const QVector3D &center, const double radius, const double highlightValue);

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1);
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v);
    virtual bool isThreadSafe();
    virtual QString getMetricNameToken();
private:
    int m_sampleRatePerKnotInterval;

    QVector<QPair<QVector3D, QVector3D> > m_boundingBoxes;
    QVector<double> m_boundingBoxHighlightValues;
    QVector<QPair<QVector3D, double> > m_boundingSpheres;
    QVector<double> m_boundingSphereHighlightValues;
};

#endif // SUBDIVISIONMANUALHIGHLIGHTMETRIC_H

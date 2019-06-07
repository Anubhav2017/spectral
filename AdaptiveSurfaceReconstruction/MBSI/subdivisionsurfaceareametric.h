#ifndef SUBDIVISIONSURFACEAREAMETRIC_H
#define SUBDIVISIONSURFACEAREAMETRIC_H

#include "subdivisionabstractmetric.h"

class SubdivisionSurfaceAreaMetric : public SubdivisionAbstractMetric
{
public:
    SubdivisionSurfaceAreaMetric(const int sampleRatePerKnotInterval);

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1);
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v);
    virtual bool isThreadSafe();
    virtual QString getMetricNameToken();
private:
    int m_sampleRatePerKnotInterval;
};

#endif // SUBDIVISIONSURFACEAREAMETRIC_H

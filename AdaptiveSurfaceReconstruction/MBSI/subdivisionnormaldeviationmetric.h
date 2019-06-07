#ifndef SUBDIVISIONNORMALDEVIATIONMETRIC_H
#define SUBDIVISIONNORMALDEVIATIONMETRIC_H

#include "subdivisionabstractmetric.h"

class SubdivisionNormalDeviationMetric : public SubdivisionAbstractMetric
{
public:
    SubdivisionNormalDeviationMetric(const int sampleRatePerKnotInterval);

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1);
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v);
    virtual bool isThreadSafe();
    virtual QString getMetricNameToken();

private:
    int m_sampleRatePerKnotInterval;
};

#endif // SUBDIVISIONNORMALDEVIATIONMETRIC_H

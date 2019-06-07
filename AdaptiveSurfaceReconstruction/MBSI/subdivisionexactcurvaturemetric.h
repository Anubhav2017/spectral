#ifndef SUBDIVISIONEXACTCURVATUREMETRIC_H
#define SUBDIVISIONEXACTCURVATUREMETRIC_H

#include "subdivisionabstractmetric.h"

class SubdivisionExactCurvatureMetric : public SubdivisionAbstractMetric
{
public:
    SubdivisionExactCurvatureMetric(const int sampleRatePerKnotInterval, const bool integrateIn3D);

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1);
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v);
    virtual bool isThreadSafe();
    virtual QString getMetricNameToken();
private:
    int m_sampleRatePerKnotInterval;
    bool m_integrateIn3D;
};

#endif // SUBDIVISIONEXACTCURVATUREMETRIC_H

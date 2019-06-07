#ifndef SUBDIVISIONABSTRACTMETRIC_H
#define SUBDIVISIONABSTRACTMETRIC_H

#include "BSpline/bsplinesurface.h"

#include <QString>

class SubdivisionAbstractMetric
{
public:
    SubdivisionAbstractMetric()
    {}

    virtual ~SubdivisionAbstractMetric()
    {}

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1) = 0;
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v) = 0;
    virtual bool isThreadSafe() = 0;
    virtual QString getMetricNameToken() = 0;
};

#endif // SUBDIVISIONABSTRACTMETRIC_H

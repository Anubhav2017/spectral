#ifndef SURFACESUBDIVIDER_H
#define SURFACESUBDIVIDER_H

#include "BSpline/bsplinesurface.h"
#include "MBSI/subdivisionabstractmetric.h"
#include "MBSI/subsurfacetree.h"

#include <QVector>

class SurfaceSubdivider
{
public:
    //Subdivision methods
    static QVector<SubSurfaceTree *> subdivideSurfacesThresholdBased(const QVector<BSplineSurface *> &surfaces, SubdivisionAbstractMetric *metric, const double threshold, const int maxDepth);
    static QVector<SubSurfaceTree *> subdivideSurfacesThresholdBased(const QVector<BSplineSurface *> &surfaces, const QVector<SubdivisionAbstractMetric *> &metrics, const QVector<double> &thresholds, const int maxDepth);

    static QVector<SubSurfaceTree *> subdivideSurfacesUntilFixedNumber(const QVector<BSplineSurface *> &surfaces, SubdivisionAbstractMetric *metric, const int targetNumberOfLeaves, const int maxDepth);
    static QVector<SubSurfaceTree *> subdivideSurfacesUntilFixedNumber(const QVector<BSplineSurface *> &surfaces, const QVector<SubdivisionAbstractMetric *> &metrics, const QVector<double> &weights, const int targetNumberOfLeaves, const int maxDepth);
};

#endif // SURFACESUBDIVIDER_H

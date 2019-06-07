#ifndef PARAMETEROPTIMIZER_H
#define PARAMETEROPTIMIZER_H

#include <QVector>
#include <QVector2D>
#include <QVector3D>

#include "BSpline/bsplinepatchnetwork.h"

class ParameterOptimizer
{
public:
    static double optimizeParameterizationFullSurfaceNetwork(BSplinePatchNetwork *bSplineSurfacePatches,
                                                             const bool alsoOptimizeBoundaries,
                                                             const bool useNonlinearSolver = true);

    static void optimizeParameterizationFullCurveNetwork(BSplinePatchNetwork *bSplineSurfaceNetwork,
                                                         const bool useNonlinearSolver = true);

    static QVector<double> optimizeParameterizationSingleCurve(BSplineCurve * const curve,
                                                               const QVector<QVector3D> &dataPoints,
                                                               const QVector<double> &parameters,
                                                               const bool keepEndpointsFixed = true,
                                                               const bool useNonlinearSolver = true);

    static QVector2D optimizeParameter_surfacePoint_nlo_bounded(BSplineSurface *surface,
                                                                const QVector3D &point,
                                                                const QVector2D &initialParameter = QVector2D(0.5, 0.5),
                                                                const int nloptMaxIterations = 10000,
                                                                const double nloptPrecision = 0.001);
    static QVector2D optimizeParameter_surfacePoint_projection_unbounded(BSplineSurface *surface,
                                                                         const QVector3D &point,
                                                                         const QVector2D &initialParameter = QVector2D(0.5, 0.5));
    static QVector2D optimizeParameter_surfacePoint_projection_bounded(BSplineSurface *surface,
                                                                       const QVector3D &point,
                                                                       const QVector2D &initialParameter = QVector2D(0.5, 0.5));

    static double optimizeParameter_curvePoint_nlo_bounded(BSplineCurve *curve,
                                                           const QVector3D &point,
                                                           const double initialParameter = 0.5,
                                                           const int nloptMaxIterations = 10000,
                                                           const double nloptPrecision = 0.001);
    static double optimizeParameter_curvePoint_projection_unbounded(BSplineCurve *curve,
                                                                    const QVector3D &point,
                                                                    const double initialParameter = 0.5);
    static double optimizeParameter_curvePoint_projection_bounded(BSplineCurve *curve,
                                                                  const QVector3D &point,
                                                                  const double initialParameter = 0.5);
};

#endif // PARAMETEROPTIMIZER_H

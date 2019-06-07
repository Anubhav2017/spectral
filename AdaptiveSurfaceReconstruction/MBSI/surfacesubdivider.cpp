#include "surfacesubdivider.h"

#include <omp.h>

QVector<SubSurfaceTree *> SurfaceSubdivider::subdivideSurfacesThresholdBased(const QVector<BSplineSurface *> &surfaces, SubdivisionAbstractMetric *metric, const double threshold, const int maxDepth)
{
    QVector<SubdivisionAbstractMetric *> metrics(1);
    QVector<double> thresholds(1);
    metrics[0] = metric;
    thresholds[0] = threshold;

    return SurfaceSubdivider::subdivideSurfacesThresholdBased(surfaces, metrics, thresholds, maxDepth);
}

QVector<SubSurfaceTree *> SurfaceSubdivider::subdivideSurfacesThresholdBased(const QVector<BSplineSurface *> &surfaces, const QVector<SubdivisionAbstractMetric *> &metrics, const QVector<double> &thresholds, const int maxDepth)
{
    const int numberOfSurfaces = surfaces.size();

    //check if one metric is not thread safe
    bool metricsAreThreadSafe = true;
    foreach (SubdivisionAbstractMetric * metric, metrics)
        if (!metric->isThreadSafe())
            metricsAreThreadSafe = false;

    const int numberOfThreads = omp_get_max_threads();
    if (!metricsAreThreadSafe)
        omp_set_num_threads(1);

    //Compute the subdivision tree
    QVector<SubSurfaceTree *> subdivisions(numberOfSurfaces);
    #pragma omp parallel for
    for (int i = 0; i < numberOfSurfaces; i++) {
        subdivisions[i] = new SubSurfaceTree(surfaces[i], metrics);
        subdivisions[i]->subdivideCentralRecursively(thresholds, maxDepth);
    }

    //If metric was not thread safe, restore original multithreading
    if (!metricsAreThreadSafe)
        omp_set_num_threads(numberOfThreads);

    return subdivisions;
}

QVector<SubSurfaceTree *> SurfaceSubdivider::subdivideSurfacesUntilFixedNumber(const QVector<BSplineSurface *> &surfaces, SubdivisionAbstractMetric *metric, const int targetNumberOfLeaves, const int maxDepth)
{
    QVector<SubdivisionAbstractMetric *> metrics(1);
    QVector<double> weights(1);
    metrics[0] = metric;
    weights[0] = 1.0;

    return SurfaceSubdivider::subdivideSurfacesUntilFixedNumber(surfaces, metrics, weights, targetNumberOfLeaves, maxDepth);
}

QVector<SubSurfaceTree *> SurfaceSubdivider::subdivideSurfacesUntilFixedNumber(const QVector<BSplineSurface *> &surfaces, const QVector<SubdivisionAbstractMetric *> &metrics, const QVector<double> &weights, const int targetNumberOfLeaves, const int maxDepth)
{
    const int numberOfSurfaces = surfaces.size();

    //check if one metric is not thread safe
    bool metricsAreThreadSafe = true;
    foreach (SubdivisionAbstractMetric * metric, metrics)
        if (!metric->isThreadSafe())
            metricsAreThreadSafe = false;

    const int numberOfThreads = omp_get_max_threads();
    if (!metricsAreThreadSafe)
        omp_set_num_threads(1);

    //Set up surfaces without subdivision
    QVector<SubSurfaceTree *> subdivisions(numberOfSurfaces);
    #pragma omp parallel for
    for (int i = 0; i < numberOfSurfaces; i++) {
        subdivisions[i] = new SubSurfaceTree(surfaces[i], metrics);
    }

    //Subdivide surface with highest value until target amount of leaves is reached
    QList<SubSurfaceTree *> leaves = subdivisions.toList();
    int numberOfLeaves = leaves.size();

    while (numberOfLeaves < targetNumberOfLeaves) {
        //find leaf with highest value
        int highestValueId = -1;
        double highestValue = -99999;
        for (int id = 0; id < numberOfLeaves; id++) {
            const double value = leaves[id]->getMetricValueSum(weights);
            if (value > highestValue && leaves[id]->getDepth() < maxDepth) {
                highestValue = value;
                highestValueId = id;
            }
        }

        SubSurfaceTree *highestValueLeaf = leaves[highestValueId];
        leaves.removeAt(highestValueId);

        highestValueLeaf->subdivideCentral();
        leaves.append(highestValueLeaf->getChildrenP()->toList());

        numberOfLeaves = leaves.size();
    }

    //If metric was not thread safe, restore original multithreading
    if (!metricsAreThreadSafe)
        omp_set_num_threads(numberOfThreads);

    return subdivisions;
}

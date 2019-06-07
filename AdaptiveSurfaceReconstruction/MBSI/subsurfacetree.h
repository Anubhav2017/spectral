#ifndef SUBSURFACETREE_H
#define SUBSURFACETREE_H

#include "subdivisionabstractmetric.h"
#include "BSpline/bsplinesurface.h"

#include <QVector>

class SubSurfaceTree
{
public:
    //Simplified constructor that defines the subsurface to span the entire base surface
    SubSurfaceTree(BSplineSurface *surf, const QVector<SubdivisionAbstractMetric *> &metrics);
    //Constructor for actual subsurface (i.e., arbitrary parts of the base surface)
    SubSurfaceTree(BSplineSurface *surf, const QVector<SubdivisionAbstractMetric *> &metrics,
               const double u0, const double u1, const double v0, const double v1, const int depth = 0, SubSurfaceTree *parent = 0);

    //Destructor (deletes all children)
    ~SubSurfaceTree();

    //Subdivision
    void subdivideCentral();
    void subdivideCentralRecursively(const QVector<double> &thresholds, const int maxDepth);

    //Metrics
    bool allMetricsBelowThresholds(const QVector<double> &thresholds) const;
    const QVector<double> &getMetricValues() const;
    double getMetricValueSum(const QVector<double> &weights) const;

    //Parameter ranges and depth
    BSplineSurface *getBaseSurface() const;
    int getDepth() const;
    double getU0() const;
    double getU1() const;
    double getUCenter() const;
    double getV0() const;
    double getV1() const;
    double getVCenter() const;

    //Tree operations
    SubSurfaceTree *getRoot();
    SubSurfaceTree *getParent() const;
    QVector<SubSurfaceTree *> getChildren() const;
    const QVector<SubSurfaceTree *> *getChildrenP() const;
    QVector<SubSurfaceTree *> getAllLeaves();
    bool isLeaf() const;

private:
    BSplineSurface *m_baseSurface;
    const QVector<SubdivisionAbstractMetric *> m_metrics;
    QVector<double> m_metricValues;

    const double m_u0;
    const double m_u1;
    const double m_v0;
    const double m_v1;
    const int m_depth;

    SubSurfaceTree *m_parent;
    QVector<SubSurfaceTree *> m_children;
};

#endif // SUBSURFACETREE_H

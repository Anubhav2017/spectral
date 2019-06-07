#include "subsurfacetree.h"

SubSurfaceTree::SubSurfaceTree(BSplineSurface *surf, const QVector<SubdivisionAbstractMetric *> &metrics):
    SubSurfaceTree(surf, metrics,
               surf->getPropertiesU()->getInnerKnotsStart(), surf->getPropertiesU()->getInnerKnotsEnd(),
               surf->getPropertiesV()->getInnerKnotsStart(), surf->getPropertiesV()->getInnerKnotsEnd(),
               0)
{

}

SubSurfaceTree::SubSurfaceTree(BSplineSurface *surf, const QVector<SubdivisionAbstractMetric *> &metrics, const double u0, const double u1, const double v0, const double v1, const int depth, SubSurfaceTree *parent):
    m_baseSurface(surf), m_metrics(metrics),
    m_u0(u0), m_u1(u1), m_v0(v0), m_v1(v1),
    m_depth(depth), m_parent(parent)
{
    const int numberOfMetrics = m_metrics.size();
    m_metricValues = QVector<double>(numberOfMetrics);

    for (int i = 0; i < numberOfMetrics; i++)
        m_metricValues[i] = m_metrics[i]->evaluate(m_baseSurface, m_u0, m_u1, m_v0, m_v1);
}

SubSurfaceTree::~SubSurfaceTree()
{
    foreach (SubSurfaceTree *surf, m_children)
        delete surf;
}

void SubSurfaceTree::subdivideCentral()
{
    const double uCenter = (m_u0 + m_u1)/2;
    const double vCenter = (m_v0 + m_v1)/2;

    m_children = QVector<SubSurfaceTree *>(4);
    m_children[0] = new SubSurfaceTree(m_baseSurface, m_metrics, m_u0, uCenter, m_v0, vCenter, m_depth + 1);
    m_children[1] = new SubSurfaceTree(m_baseSurface, m_metrics, m_u0, uCenter, vCenter, m_v1, m_depth + 1);
    m_children[2] = new SubSurfaceTree(m_baseSurface, m_metrics, uCenter, m_u1, m_v0, vCenter, m_depth + 1);
    m_children[3] = new SubSurfaceTree(m_baseSurface, m_metrics, uCenter, m_u1, vCenter, m_v1, m_depth + 1);
}

void SubSurfaceTree::subdivideCentralRecursively(const QVector<double> &thresholds, const int maxDepth)
{
    if (!this->allMetricsBelowThresholds(thresholds) && m_depth < maxDepth) {
        this->subdivideCentral();
        foreach (SubSurfaceTree *child, m_children)
            child->subdivideCentralRecursively(thresholds, maxDepth);
    }
}

bool SubSurfaceTree::allMetricsBelowThresholds(const QVector<double> &thresholds) const
{
    const int numberOfMetrics = m_metricValues.size();
    for (int i = 0; i < numberOfMetrics; i++)
        if (m_metricValues[i] > thresholds[i])
            return false;

    return true;
}

const QVector<double> &SubSurfaceTree::getMetricValues() const
{
    return m_metricValues;
}

double SubSurfaceTree::getMetricValueSum(const QVector<double> &weights) const
{
    const int numberOfMetrics = m_metricValues.size();
    double result = 0;
    for (int i = 0; i < numberOfMetrics; i++)
        result += m_metricValues[i] * weights[i];

    return result;
}

BSplineSurface *SubSurfaceTree::getBaseSurface() const
{
    return m_baseSurface;
}

int SubSurfaceTree::getDepth() const
{
    return m_depth;
}

double SubSurfaceTree::getU0() const
{
    return m_u0;
}

double SubSurfaceTree::getU1() const
{
    return m_u1;
}

double SubSurfaceTree::getUCenter() const
{
    return (m_u0 + m_u1)/2;
}

double SubSurfaceTree::getV0() const
{
    return m_v0;
}

double SubSurfaceTree::getV1() const
{
    return m_v1;
}

double SubSurfaceTree::getVCenter() const
{
    return (m_v0 + m_v1)/2;
}

SubSurfaceTree *SubSurfaceTree::getRoot()
{
    if (!m_parent || m_parent == this)
        return this;
    else
        return m_parent->getRoot();
}

SubSurfaceTree *SubSurfaceTree::getParent() const
{
    return m_parent;
}

QVector<SubSurfaceTree *> SubSurfaceTree::getChildren() const
{
    return m_children;
}

const QVector<SubSurfaceTree *> *SubSurfaceTree::getChildrenP() const
{
    return &m_children;
}

QVector<SubSurfaceTree *> SubSurfaceTree::getAllLeaves()
{
    QVector<SubSurfaceTree *> result;

    if (this->isLeaf()) {
        result << this;
    } else {
        foreach (SubSurfaceTree *child, m_children) {
            result << child->getAllLeaves();
        }
    }

    return result;
}

bool SubSurfaceTree::isLeaf() const
{
    return m_children.isEmpty();
}

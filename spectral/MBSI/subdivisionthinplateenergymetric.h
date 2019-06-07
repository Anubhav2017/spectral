#ifndef SUBDIVISIONTHINPLATEENERGYMETRIC_H
#define SUBDIVISIONTHINPLATEENERGYMETRIC_H

#include "subdivisionabstractmetric.h"
#include "BSpline/bsplinepatchnetwork.h"
#include "energymatrixgenerator.h"

#include <Sparse>

#include <QHash>
#include <QString>

class SubdivisionThinPlateEnergyMetric : public SubdivisionAbstractMetric
{
public:
    SubdivisionThinPlateEnergyMetric(EnergyMatrixGenerator *energyMatGenerator, BSplineProperties *properties);
    virtual ~SubdivisionThinPlateEnergyMetric();

    virtual double evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1);
    virtual double evaluatePoint(BSplineSurface *surf, const double u, const double v);
    virtual bool isThreadSafe();    //this metric is not thread safe (potential conflict when matrices get inserted into hash)
    virtual QString getMetricNameToken();

    void precomputeThinPlateMatrices(const int subdivisionDepth);

    Eigen::SparseMatrix<double> *getThinPlateEnergyMatrixForRegion(const double u0, const double u1, const double v0, const double v1);
private:
    void precomputeRecursive(const int depth, const int maxDepth, const double u0, const double u1, const double v0, const double v1);

    EnergyMatrixGenerator *m_energyMatrixGenerator;
    QHash<QString, Eigen::SparseMatrix<double> *> m_energyMatrices;

    BSplineProperties *m_properties;
};

#endif // SUBDIVISIONTHINPLATEENERGYMETRIC_H

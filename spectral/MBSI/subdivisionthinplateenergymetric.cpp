#include "subdivisionthinplateenergymetric.h"

#include "splineevaluator.h"

SubdivisionThinPlateEnergyMetric::SubdivisionThinPlateEnergyMetric(EnergyMatrixGenerator *energyMatGenerator, BSplineProperties *properties) :
    SubdivisionAbstractMetric(),
    m_energyMatrixGenerator(energyMatGenerator),
    m_properties(properties)
{

}

SubdivisionThinPlateEnergyMetric::~SubdivisionThinPlateEnergyMetric()
{
    foreach(Eigen::SparseMatrix<double> *mat, m_energyMatrices.values())
        delete mat;
}

double SubdivisionThinPlateEnergyMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
{
    //calculate patch energy in square from [u0, v0] to [u1, v1]
    Eigen::SparseMatrix<double> energyMatrix = *this->getThinPlateEnergyMatrixForRegion(u0, u1, v0, v1);
    return SplineEvaluator::calculatePatchEnergy(surf, energyMatrix);
}

double SubdivisionThinPlateEnergyMetric::evaluatePoint(BSplineSurface *surf, const double u, const double v)
{
    return surf->evaluateDerivative(u, v, 2, 0).lengthSquared() + (2*surf->evaluateDerivative(u, v, 1, 1).lengthSquared()) + surf->evaluateDerivative(u, v, 0, 2).lengthSquared();
}

bool SubdivisionThinPlateEnergyMetric::isThreadSafe()
{
    return false;
}

QString SubdivisionThinPlateEnergyMetric::getMetricNameToken()
{
    return QString("ThinPlate");
}

void SubdivisionThinPlateEnergyMetric::precomputeThinPlateMatrices(const int subdivisionDepth)
{
    const double intervalStart = m_properties->getInnerKnotsStart();
    const double intervalEnd = m_properties->getInnerKnotsEnd();
    this->precomputeRecursive(0, subdivisionDepth, intervalStart, intervalEnd, intervalStart, intervalEnd);
}

Eigen::SparseMatrix<double> *SubdivisionThinPlateEnergyMetric::getThinPlateEnergyMatrixForRegion(const double u0, const double u1, const double v0, const double v1)
{
    QString region = QString("%1%2%3%4").arg(u0).arg(u1).arg(v0).arg(v1);
    if (m_energyMatrices.contains(region))
        return m_energyMatrices[region];
    else {
        const int numTotalCP = m_properties->getNumberOfControlPoints() * m_properties->getNumberOfControlPoints();
        Eigen::SparseMatrix<double> *newMatrix = new Eigen::SparseMatrix<double>(numTotalCP, numTotalCP);
        (*newMatrix) += m_energyMatrixGenerator->getEnergyMatrixSurface(m_properties, 2, 0, u0, u1, v0, v1);
        (*newMatrix) += m_energyMatrixGenerator->getEnergyMatrixSurface(m_properties, 1, 1, u0, u1, v0, v1);
        (*newMatrix) += m_energyMatrixGenerator->getEnergyMatrixSurface(m_properties, 0, 2, u0, u1, v0, v1);
        m_energyMatrices[region] = newMatrix;
        return newMatrix;
    }
}

void SubdivisionThinPlateEnergyMetric::precomputeRecursive(const int depth, const int maxDepth, const double u0, const double u1, const double v0, const double v1)
{
    getThinPlateEnergyMatrixForRegion(u0, u1, v0, v1);
    if (depth < maxDepth) {
        const double uCenter = 0.5 * (u0 + u1);
        const double vCenter = 0.5 * (v0 + v1);
        this->precomputeRecursive(depth + 1, maxDepth, u0, uCenter, v0, vCenter);
        this->precomputeRecursive(depth + 1, maxDepth, u0, uCenter, vCenter, v1);
        this->precomputeRecursive(depth + 1, maxDepth, uCenter, u1, v0, vCenter);
        this->precomputeRecursive(depth + 1, maxDepth, uCenter, u1, vCenter, v1);
    }
}


#include "subdivisionexactcurvaturemetric.h"

SubdivisionExactCurvatureMetric::SubdivisionExactCurvatureMetric(const int sampleRatePerKnotInterval, const bool integrateIn3D):
    SubdivisionAbstractMetric(),
    m_sampleRatePerKnotInterval(sampleRatePerKnotInterval), m_integrateIn3D(integrateIn3D)
{

}

double SubdivisionExactCurvatureMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
{
    const double uDiff = u1 - u0;
    const double vDiff = v1 - v0;

    const int sampleRateU = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesU()->getNumberOfKnotsInOpenInterval(u0, u1));
    const int sampleRateV = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesV()->getNumberOfKnotsInOpenInterval(v0, v1));

    double result = 0;
    const double areaElementParameterSpace = uDiff * vDiff / (sampleRateU * sampleRateV); //needed for the numeric integration
    for (int iu = 0; iu < sampleRateU; iu++) {
        const double u = u0 + uDiff * (((double) iu + 0.5)/ sampleRateU);
        for (int iv = 0; iv < sampleRateV; iv++) {
            const double v = v0 + vDiff * (((double) iv + 0.5)/ sampleRateV);
            const double k1k1_k2k2 = this->evaluatePoint(surf, u, v);

            if (m_integrateIn3D) {
                const double surfaceElement = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v)).length();

                result += k1k1_k2k2 * surfaceElement;
            } else {
                result += k1k1_k2k2;
            }
        }
    }

    return result * areaElementParameterSpace;
}

double SubdivisionExactCurvatureMetric::evaluatePoint(BSplineSurface *surf, const double u, const double v)
{
    double K, H;
    surf->evaluateGaussianAndMeanCurvature(u, v, K, H);
    return ((4 * H * H) - (2 * K));  //equivalent to k1*k1 + k2*k2, since H = (k1+k2)/2, K = k1*k2
}

bool SubdivisionExactCurvatureMetric::isThreadSafe()
{
    return true;
}

QString SubdivisionExactCurvatureMetric::getMetricNameToken()
{
    return QString("ExactCurv");
}


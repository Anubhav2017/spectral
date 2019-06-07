#include "subdivisionsurfaceareametric.h"

SubdivisionSurfaceAreaMetric::SubdivisionSurfaceAreaMetric(const int sampleRatePerKnotInterval):
    SubdivisionAbstractMetric(),
    m_sampleRatePerKnotInterval(sampleRatePerKnotInterval)
{

}

double SubdivisionSurfaceAreaMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
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

            const double surfaceElement = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v)).length();

            result += surfaceElement;
        }
    }

    return result * areaElementParameterSpace;
}

double SubdivisionSurfaceAreaMetric::evaluatePoint(BSplineSurface *, const double, const double)
{
    return 1;
}

bool SubdivisionSurfaceAreaMetric::isThreadSafe()
{
    return true;
}

QString SubdivisionSurfaceAreaMetric::getMetricNameToken()
{
    return QString("Area");
}


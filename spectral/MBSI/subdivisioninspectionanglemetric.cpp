#include "subdivisioninspectionanglemetric.h"

#include <cmath>

SubdivisionInspectionAngleMetric::SubdivisionInspectionAngleMetric(const double cameraDistance, const int sampleRatePerKnotInterval):
    SubdivisionAbstractMetric(),
    m_sampleRatePerKnotInterval(sampleRatePerKnotInterval),
    m_cameraDistance(cameraDistance)
{

}

double SubdivisionInspectionAngleMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
{
    const double uDiff = u1 - u0;
    const double vDiff = v1 - v0;
    const int sampleRateU = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesU()->getNumberOfKnotsInOpenInterval(u0, u1));
    const int sampleRateV = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesV()->getNumberOfKnotsInOpenInterval(v0, v1));

    const double uCenter = (u0 + u1)/2;
    const double vCenter = (v0 + v1)/2;
    const QVector3D pCenter = surf->evaluate(uCenter, vCenter);
    const QVector3D normalCenter = QVector3D::crossProduct(surf->evaluateTu(uCenter, vCenter), surf->evaluateTv(uCenter, vCenter)).normalized();
    const QVector3D camPosition = pCenter + m_cameraDistance * normalCenter;

    double maxDeviation = 0;
    for (int iu = 0; iu < sampleRateU; iu++) {
        const double u = u0 + uDiff * (((double) iu + 0.5)/ sampleRateU);
        for (int iv = 0; iv < sampleRateV; iv++) {
            const double v = v0 + vDiff * (((double) iv + 0.5)/ sampleRateV);

            const QVector3D pos = surf->evaluate(u, v);
            const QVector3D normal = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v)).normalized();
            const QVector3D posToCam = (camPosition - pos).normalized();

            const double deviation = (normal - posToCam).length();
            if (deviation > maxDeviation)
                maxDeviation = deviation;
        }
    }

    return acos((2 - maxDeviation * maxDeviation)/2) * 180/M_PI;
}

double SubdivisionInspectionAngleMetric::evaluatePoint(BSplineSurface *surf, const double u, const double v)
{
    const double uCenter = (surf->getPropertiesU()->getInnerKnotsStart() + surf->getPropertiesU()->getInnerKnotsEnd())/2;
    const double vCenter = (surf->getPropertiesV()->getInnerKnotsStart() + surf->getPropertiesV()->getInnerKnotsEnd())/2;
    const QVector3D pCenter = surf->evaluate(uCenter, vCenter);
    const QVector3D normalCenter = QVector3D::crossProduct(surf->evaluateTu(uCenter, vCenter), surf->evaluateTv(uCenter, vCenter)).normalized();
    const QVector3D camPosition = pCenter + m_cameraDistance * normalCenter;

    const QVector3D pos = surf->evaluate(u, v);
    const QVector3D normal = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v)).normalized();
    const QVector3D posToCam = (camPosition - pos).normalized();

    const double deviation = (normal - posToCam).length();
    return acos((2 - deviation * deviation)/2) * 180/M_PI;
}

bool SubdivisionInspectionAngleMetric::isThreadSafe()
{
    return true;
}

QString SubdivisionInspectionAngleMetric::getMetricNameToken()
{
    return QString("InspAngle");
}

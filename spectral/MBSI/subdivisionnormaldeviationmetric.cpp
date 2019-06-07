#include "subdivisionnormaldeviationmetric.h"

#include "external/Miniball.hpp"
#include <cmath>

SubdivisionNormalDeviationMetric::SubdivisionNormalDeviationMetric(const int sampleRatePerKnotInterval):
    SubdivisionAbstractMetric(),
    m_sampleRatePerKnotInterval(sampleRatePerKnotInterval)
{

}

double SubdivisionNormalDeviationMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
{
    const double uDiff = u1 - u0;
    const double vDiff = v1 - v0;
    const int sampleRateU = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesU()->getNumberOfKnotsInOpenInterval(u0, u1));
    const int sampleRateV = m_sampleRatePerKnotInterval * (1 + surf->getPropertiesV()->getNumberOfKnotsInOpenInterval(v0, v1));

    QVector<QVector<QVector3D> > normalField(sampleRateU, QVector<QVector3D>(sampleRateV));
    QVector3D averageNormal(0, 0, 0);
    double totalArea = 0;
    for (int iu = 0; iu < sampleRateU; iu++) {
        const double u = u0 + uDiff * (((double) iu + 0.5)/ sampleRateU);
        for (int iv = 0; iv < sampleRateV; iv++) {
            const double v = v0 + vDiff * (((double) iv + 0.5)/ sampleRateV);
            QVector3D normal = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v));
            const double surfaceElement = normal.length();

            normal /= surfaceElement;

            totalArea += surfaceElement;
            averageNormal += normal * surfaceElement;

            normalField[iu][iv] = normal;
        }
    }

    averageNormal /= totalArea; //Do we need to do that? Normalization should take care of it anyway
    averageNormal.normalize();

    double maxAngle = 0;
    for (int iu = 0; iu < sampleRateU; iu++) {
        for (int iv = 0; iv < sampleRateV; iv++) {
            const double angle = acos(QVector3D::dotProduct(normalField[iu][iv], averageNormal));
            if (angle > maxAngle)
                maxAngle = angle;
        }
    }

    return maxAngle * 180/M_PI;
}

double SubdivisionNormalDeviationMetric::evaluatePoint(BSplineSurface *, const double, const double)
{
    return 0;
}

bool SubdivisionNormalDeviationMetric::isThreadSafe()
{
    return true;
}

QString SubdivisionNormalDeviationMetric::getMetricNameToken()
{
    return QString("normalDev");
}

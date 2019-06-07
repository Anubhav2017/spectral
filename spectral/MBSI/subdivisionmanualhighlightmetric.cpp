#include "subdivisionmanualhighlightmetric.h"

SubdivisionManualHighlightMetric::SubdivisionManualHighlightMetric(const int sampleRatePerKnotInterval):
    SubdivisionAbstractMetric(),
    m_sampleRatePerKnotInterval(sampleRatePerKnotInterval)
{

}

void SubdivisionManualHighlightMetric::addHighlightBoundingBox(const QVector3D &corner0, const QVector3D &corner1, const double highlightValue)
{
    float lowerX = corner0.x();
    float upperX = corner1.x();
    if (lowerX > upperX) {
        lowerX = corner1.x();
        upperX = corner0.x();
    }

    float lowerY = corner0.y();
    float upperY = corner1.y();
    if (lowerY > upperY) {
        lowerY = corner1.y();
        upperY = corner0.y();
    }

    float lowerZ = corner0.z();
    float upperZ = corner1.z();
    if (lowerZ > upperZ) {
        lowerZ = corner1.z();
        upperZ = corner0.z();
    }

    m_boundingBoxes.append(QPair<QVector3D, QVector3D>(QVector3D(lowerX, lowerY, lowerZ), QVector3D(upperX, upperY, upperZ)));
    m_boundingBoxHighlightValues.append(highlightValue);
}

void SubdivisionManualHighlightMetric::addHighlightBoundingSphere(const QVector3D &center, const double radius, const double highlightValue)
{
    m_boundingSpheres.append(QPair<QVector3D, double>(QVector3D(center), radius));
    m_boundingSphereHighlightValues.append(highlightValue);
}

double SubdivisionManualHighlightMetric::evaluate(BSplineSurface *surf, const double u0, const double u1, const double v0, const double v1)
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
            const double highlightValue = this->evaluatePoint(surf, u, v);

            if (highlightValue != 0) {
                const double surfaceElement = QVector3D::crossProduct(surf->evaluateTu(u, v), surf->evaluateTv(u, v)).length();

                result += surfaceElement * highlightValue;
            }
        }
    }

    return result * areaElementParameterSpace;
}

double SubdivisionManualHighlightMetric::evaluatePoint(BSplineSurface *surf, const double u, const double v)
{
    const QVector3D point = surf->evaluate(u, v);

    const int numBoundingBoxes = m_boundingBoxes.size();
    double result = 0;
    for (int i = 0; i < numBoundingBoxes; i++) {
        const QVector3D &corner0 = m_boundingBoxes[i].first;
        const QVector3D &corner1 = m_boundingBoxes[i].second;

        if (point.x() >= corner0.x() && point.x() <= corner1.x()
            && point.y() >= corner0.y() && point.y() <= corner1.y()
            && point.z() >= corner0.z() && point.z() <= corner1.z())
            result += m_boundingBoxHighlightValues[i];
    }

    const int numBoundingSpheres = m_boundingSpheres.size();
    for (int i = 0; i < numBoundingSpheres; i++) {
        const QVector3D &center = m_boundingSpheres[i].first;
        const double rSquared = m_boundingSpheres[i].second * m_boundingSpheres[i].second;
        //squaring the equation, gets rid of the computation of the squareroot

        if ((point - center).lengthSquared() < rSquared)
            result += m_boundingSphereHighlightValues[i];
    }

    return result;
}

bool SubdivisionManualHighlightMetric::isThreadSafe()
{
    return true;
}

QString SubdivisionManualHighlightMetric::getMetricNameToken()
{
    return QString("ManHighlight");
}


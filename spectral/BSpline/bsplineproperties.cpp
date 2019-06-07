#include "bsplineproperties.h"

BSplineProperties *BSplineProperties::createInstanceUniformKnotVectorMultipleEnds(const int order, const int numberOfControlPoints, const double tStart, const double tEnd)
{
    const double tDiff = tEnd - tStart;
    //At least order+1 control points required
    if (numberOfControlPoints < order + 1)
        return 0;

    const int numberOfKnots = numberOfControlPoints + order + 1;
    QVector<double> knots(numberOfKnots);

    //First order + 1 knots are tStart
    for (int i = 0; i < order + 1; i++)
        knots[i] = tStart;

    //Inner knots are equally distributed
    const int numberOfInnerKnots = numberOfControlPoints - order - 1;   //numCtrPoints + order + 1 - 2*(order +1)
    for (int i = 0; i < numberOfInnerKnots; i++)
        knots[i + order + 1] = tStart + (double) (i+1) * tDiff / (numberOfInnerKnots + 1);

    //Last order + 1 knots are tEnd
    for (int i = order + 1 + numberOfInnerKnots; i < numberOfKnots; i++)
        knots[i] = tEnd;

    return new BSplineProperties(order, numberOfControlPoints, knots);
}

double BSplineProperties::N(const double t, const int i) const
{
    return N(t, i, m_order);
}

double BSplineProperties::N(const double t, const int i, const int k) const
{
    if (t < m_knots[i] || (i + k + 1 < m_knots.size() && t > m_knots[i + k + 1])) {
        return 0;
    }

    //This part is not nice, but helps to avoid the problem that N is only defined on the open interval [start, end[
    if (t == m_knots[i + 1] && t == m_knots.last()) {
        return 1;
    }

    if (k == 0) {
        if (m_knots[i] <= t && t < m_knots[i + 1])
            return 1;
        else
            return 0;
    } else {
        if (m_knots[i] == m_knots[i + k + 1])
            return 1;

        const double N1 = N(t, i  , k-1);
        const double N2 = N(t, i+1, k-1);
        double firstHalf = 0;
        double secondHalf = 0;

        if ((m_knots[i + k] - m_knots[i]) != 0)
            firstHalf = N1 * (t - m_knots[i]) / (m_knots[i + k] - m_knots[i]);
        if (((m_knots[i + k + 1] - m_knots[i + 1])) != 0)
            secondHalf = N2 * (m_knots[i + k + 1] - t) / (m_knots[i + k + 1] - m_knots[i + 1]);

        return firstHalf + secondHalf;
    }
}

double BSplineProperties::NDerivative(const double t, const int i, const int derivative) const
{
    return this->NDerivative(t, i, derivative, m_order);
}

double BSplineProperties::NDerivative(const double t, const int i, const int derivative, const int k) const
{
    const double firstInterval = m_knots[i + k] - m_knots[i];
    const double secondInterval = m_knots[i + k + 1] - m_knots[i + 1];
    double firstHalf = 0;
    double secondHalf = 0;
    if (derivative == 0) {
        return N(t, i, k);
    } else if (derivative == 1) {
        if (firstInterval != 0)
            firstHalf = N(t, i, k - 1) * k / firstInterval;
        if (secondInterval != 0)
            secondHalf = N(t, i+1, k - 1) * k / secondInterval;
    } else if (derivative <= k) {   //The first derivative of k = 0 (constant functions) would be 0 -> derivative > k equals 0
        if (firstInterval != 0)
            firstHalf = NDerivative(t, i, derivative - 1, k - 1) * k / firstInterval;
        if (secondInterval != 0)
            secondHalf = NDerivative(t, i+1, derivative - 1, k - 1) * k / secondInterval;
    } else {
        return 0;
    }

    return firstHalf - secondHalf;
}

bool BSplineProperties::equalTo(const BSplineProperties &other) const
{
    if (m_order != other.getOrder() || m_numberOfControlPoints != other.getNumberOfControlPoints())
        return false;

    const int numKnots = m_knots.size();
    if (other.getKnots().size() != numKnots)
        return false;

    for (int i = 0; i < numKnots; i++)
        if (other.getKnotsP()->at(i) != m_knots[i])
            return false;

    return true;
}

int BSplineProperties::getNumberOfKnotsInOpenInterval(const double tStart, const double tEnd) const
{
    int result = 0;
    foreach (const double knot, m_knots)
        if (knot > tStart && knot < tEnd)
            result++;

    return result;
}

int BSplineProperties::getOrder() const
{
    return m_order;
}

int BSplineProperties::getNumberOfControlPoints() const
{
    return m_numberOfControlPoints;
}

QVector<double> BSplineProperties::getKnots() const
{
    return m_knots;
}

const QVector<double> *BSplineProperties::getKnotsP() const
{
    return &m_knots;
}

QVector<double> BSplineProperties::getInnerKnotVector() const
{
    const int start = this->getFirstInnerKnotIndex();
    const int end = this->getLastInnerKnotIndex();
    QVector<double> innerKnots(1 + end - start);
    for (int i = start; i <= end; i++)
        innerKnots[i - start] = this->getKnotsP()->at(i);

    return innerKnots;
}

int BSplineProperties::getFirstInnerKnotIndex() const
{
    return m_order;
}

int BSplineProperties::getLastInnerKnotIndex() const
{
    return m_numberOfControlPoints;
}

double BSplineProperties::getInnerKnotsStart() const
{
    return m_knots[this->getFirstInnerKnotIndex()];
}

double BSplineProperties::getInnerKnotsEnd() const
{
    return m_knots[this->getLastInnerKnotIndex()];
}

BSplineProperties::BSplineProperties(const int order, const int numberOfControlPoints, const QVector<double> &knots):
    m_order(order), m_numberOfControlPoints(numberOfControlPoints), m_knots(knots)
{

}

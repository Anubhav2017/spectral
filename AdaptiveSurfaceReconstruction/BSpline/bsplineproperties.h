#ifndef BSPLINEPROPERTIES_H
#define BSPLINEPROPERTIES_H

#include <QVector>

class BSplineProperties
{
public:
    //Setup methods
    static BSplineProperties *createInstanceUniformKnotVectorMultipleEnds(const int order, const int numberOfControlPoints, const double tStart, const double tEnd);

    //Basis functions
    double N(const double t, const int i) const;
    double N(const double t, const int i, const int k) const;
    double NDerivative(const double t, const int i, const int derivative) const;
    double NDerivative(const double t, const int i, const int derivative, const int k) const;

    //Other functions
    bool equalTo(const BSplineProperties &other) const;
    int getNumberOfKnotsInOpenInterval(const double tStart, const double tEnd) const;

    //Access functions
    int getOrder() const;
    int getNumberOfControlPoints() const;
    QVector<double> getKnots() const;
    const QVector<double> *getKnotsP() const;
    QVector<double> getInnerKnotVector() const;
    int getFirstInnerKnotIndex() const;
    int getLastInnerKnotIndex() const;
    double getInnerKnotsStart() const;
    double getInnerKnotsEnd() const;
private:
    //Constructor is only supposed to be used by static factory methods
    BSplineProperties(const int order, const int numberOfControlPoints, const QVector<double> &knots);

    //Defining properties of the B-spline
    const int m_order;
    const int m_numberOfControlPoints;
    const QVector<double> m_knots;
};

#endif // BSPLINEPROPERTIES_H

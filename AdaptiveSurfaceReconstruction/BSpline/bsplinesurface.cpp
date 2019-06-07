#include "bsplinesurface.h"

#include <cmath>

BSplineSurface::BSplineSurface(const BSplineProperties * const propertiesU, const BSplineProperties * const propertiesV, const QVector<QVector<QVector3D> > &controlPoints):
    m_propertiesU(propertiesU), m_propertiesV(propertiesV), m_controlPoints(controlPoints)
{

}

QVector3D BSplineSurface::evaluate(double u, double v)
{
    /*TODO
    / for large data, B=0 in many cases
    / -> only consider those i,j for which x in [knot_i - ?, knot_i + ?] x [...]
    */
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++)
        for (int j = 0; j < m_controlPoints[i].size(); j++)
            result += m_controlPoints[i][j] * m_propertiesV->N(v, j) * m_propertiesU->N(u, i);

    return result;
}

QVector3D BSplineSurface::evaluateTu(double u, double v)
{
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++)
        for (int j = 0; j < m_controlPoints[i].size(); j++)
            result += m_controlPoints[i][j] * m_propertiesV->N(v, j) * m_propertiesU->NDerivative(u, i, 1);

    return result;
}

QVector3D BSplineSurface::evaluateTv(double u, double v)
{
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++)
        for (int j = 0; j < m_controlPoints[i].size(); j++)
            result += m_controlPoints[i][j] * m_propertiesV->NDerivative(v, j, 1) * m_propertiesU->N(u, i);

    return result;
}

QVector3D BSplineSurface::evaluateDerivative(double u, double v, int du, int dv)
{
    QVector3D result(0, 0, 0);
    for (int i = 0; i < m_controlPoints.size(); i++)
        for (int j = 0; j < m_controlPoints[i].size(); j++)
            result += m_controlPoints[i][j] * m_propertiesV->NDerivative(v, j, dv) * m_propertiesU->NDerivative(u, i, du);

    return result;
}

void BSplineSurface::evaluateFirstFundamentalForm(double u, double v, double &E, double &F, double &G)
{
    const QVector3D Tu = this->evaluateTu(u, v);
    const QVector3D Tv = this->evaluateTv(u, v);

    E = QVector3D::dotProduct(Tu, Tu);
    F = QVector3D::dotProduct(Tu, Tv);
    G = QVector3D::dotProduct(Tv, Tv);
}

void BSplineSurface::evaluateSecondFundamentalForm(double u, double v, double &L, double &M, double &N)
{
    const QVector3D normal = QVector3D::crossProduct(this->evaluateTu(u, v), this->evaluateTv(u, v)).normalized();
    const QVector3D Tuu = this->evaluateDerivative(u, v, 2, 0);
    const QVector3D Tuv = this->evaluateDerivative(u, v, 1, 1);
    const QVector3D Tvv = this->evaluateDerivative(u, v, 0, 2);

    L = QVector3D::dotProduct(Tuu, normal);
    M = QVector3D::dotProduct(Tuv, normal);
    N = QVector3D::dotProduct(Tvv, normal);
}

double BSplineSurface::evaluateGaussianCurvature(double u, double v)
{
    double E, F, G; //First fundamental form
    double L, M, N; //Second fundamental form
    this->evaluateFirstFundamentalForm(u, v, E, F, G);
    this->evaluateSecondFundamentalForm(u, v, L, M, N);

    if (E*G == F*F) //degenerate case, happens when Tu == -Tv -> division by 0 not really defined
        return 0;

    return (L*N - M*M)/(E*G - F*F);
}

double BSplineSurface::evaluateMeanCurvature(double u, double v)
{
    double E, F, G; //First fundamental form
    double L, M, N; //Second fundamental form
    this->evaluateFirstFundamentalForm(u, v, E, F, G);
    this->evaluateSecondFundamentalForm(u, v, L, M, N);

    if (E*G == F*F) //degenerate case, happens when Tu == -Tv -> division by 0 not really defined
        return 0;

    return (E*N - 2*F*M + G*L)/(2*(E*G - F*F));
}

void BSplineSurface::evaluateGaussianAndMeanCurvature(double u, double v, double &K, double &H)
{
    double E, F, G; //First fundamental form
    double L, M, N; //Second fundamental form
    this->evaluateFirstFundamentalForm(u, v, E, F, G);
    this->evaluateSecondFundamentalForm(u, v, L, M, N);

    if (E*G - F*F == 0) { //degenerate case, happens when Tu == -Tv -> division by 0 not really defined
        K = 0;
        H = 0;
    } else {
        K = (L*N - M*M)/(E*G - F*F);
        H = (E*N - 2*F*M + G*L)/(2*(E*G - F*F));
    }
}

void BSplineSurface::evaluatePrincipalCurvatures(double u, double v, double &k1, double &k2)
{
    double K, H;
    this->evaluateGaussianAndMeanCurvature(u, v, K, H);

    const double sqrtPart = sqrt(H*H - K);

    k1 = H + sqrtPart;
    k2 = H - sqrtPart;
}

void BSplineSurface::translate(const QVector3D translationVector)
{
    const int numCtrPoints = m_controlPoints.size();
    for (int i = 0; i < numCtrPoints; i++)
        for (int j = 0; j < numCtrPoints; j++)
            m_controlPoints[i][j] += translationVector;
}

void BSplineSurface::scale(double factorX, double factorY, double factorZ)
{
    const int numCtrPoints = m_controlPoints.size();
    for (int i = 0; i < numCtrPoints; i++)
        for (int j = 0; j < numCtrPoints; j++) {
            QVector3D pos = m_controlPoints[i][j];
            pos.setX(pos.x() * factorX);
            pos.setY(pos.y() * factorY);
            pos.setZ(pos.z() * factorZ);
            m_controlPoints[i][j] = pos;
        }
}

const BSplineProperties *BSplineSurface::getPropertiesU() const
{
    return m_propertiesU;
}

const BSplineProperties *BSplineSurface::getPropertiesV() const
{
    return m_propertiesV;
}

QVector<QVector<QVector3D> > *BSplineSurface::getControlPointsP()
{
    return &m_controlPoints;
}

QVector<QVector3D> BSplineSurface::getControlPointsAsVector()
{
    const int numControlPointsU = m_propertiesU->getNumberOfControlPoints();
    const int numControlPointsV = m_propertiesV->getNumberOfControlPoints();

    QVector<QVector3D> result(numControlPointsU * numControlPointsV);
    for (int i = 0; i < numControlPointsU; i++) {
        for (int j = 0; j < numControlPointsV; j++) {
            const int index = BSplineSurface::matToVecIndexLocal(i, j, numControlPointsU);
            result[index] = m_controlPoints[i][j];
        }
    }

    return result;
}

int BSplineSurface::matToVecIndexLocal(const int &i, const int &j, const int &numberOfControlPointsInUDirection)
{
    return i * numberOfControlPointsInUDirection + j;
}

int BSplineSurface::matToVecIndexLocal(const QPair<int, int> &ij, const int &numberOfControlPointsInUDirection)
{
    return ij.first * numberOfControlPointsInUDirection + ij.second;
}

void BSplineSurface::vecToMatIndexLocal(const int &iVec, int &iOut, int &jOut, const int &numberOfControlPointsInUDirection)
{
    iOut = iVec / numberOfControlPointsInUDirection;
    jOut = iVec % numberOfControlPointsInUDirection;
}

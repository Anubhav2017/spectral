#include "energymatrixgenerator.h"

#include "BSpline/bsplinesurface.h"

#include <gsl/gsl_integration.h>

#include <QStringList>
#include <QDir>
#include <QFile>
#include <QByteArray>

#include <QDebug>

EnergyMatrixGenerator::EnergyMatrixGenerator():
    m_defaultFolder(QString(".")), m_loadFromFile(false), m_saveToFile(false)
{

}

EnergyMatrixGenerator::EnergyMatrixGenerator(const QString &defaultFolder, const bool loadFromFileIfPossible, const bool saveToFileIfPossible):
    m_defaultFolder(defaultFolder), m_loadFromFile(loadFromFileIfPossible), m_saveToFile(saveToFileIfPossible)
{

}

Eigen::SparseMatrix<double> EnergyMatrixGenerator::getEnergyMatrixCurve(BSplineProperties *properties, const int derivative, const double t0, const double t1)
{
    const QVector<double> *knots = properties->getKnotsP();
    const int lastInnerKnotIndex = properties->getLastInnerKnotIndex();

    const int numberOfGaussPoints = properties->getOrder() - derivative + 5;   //only +1 needed, but we want to prevent numerical errors
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(numberOfGaussPoints);

    const int numberOfControlPoints = properties->getNumberOfControlPoints();
    Eigen::SparseMatrix<double> matrix(numberOfControlPoints, numberOfControlPoints);

    for (int i = 0; i < numberOfControlPoints; i++) {
        for (int j = 0; j < numberOfControlPoints; j++) {
            double integral = 0;
            //for (int iKnots = 0; iKnots < sizeOfKnotVector - 1; iKnots++) {
            for (int iKnots = properties->getFirstInnerKnotIndex(); iKnots < lastInnerKnotIndex; iKnots++) { //only consider inner knot vectors
                const double a = knots->at(iKnots);
                const double b = knots->at(iKnots + 1);
                if (a != b) {
                    double ta = a;
                    double tb = b;
                    if (ta < t0)
                        ta = t0;
                    if (tb > t1)
                        tb = t1;
                    if (tb > ta) {
                        for (int k = 0; k < numberOfGaussPoints; k++) {
                            double xi, wi;
                            gsl_integration_glfixed_point(ta, tb, k, &xi, &wi, table);
                            integral += wi * properties->NDerivative(xi, i, derivative) * properties->NDerivative(xi, j, derivative);
                        }
                    }
                }
            }

            matrix.insert(i, j) = integral;
        }
    }

    gsl_integration_glfixed_table_free(table);

    return matrix;
}

Eigen::SparseMatrix<double> EnergyMatrixGenerator::getEnergyMatrixSurface(BSplineProperties *properties, const int derivativeU, const int derivativeV, const double u0, const double u1, const double v0, const double v1)
{
//    const double paramA = paramInterval.first;
//    const double paramB = paramInterval.second;
//    const QString filename = m_defaultFolder + QString("/energymatrix_n%1_o%2_du%3_dv%4_in_%5_%6_%7_%8_%9_%10.mat.dm").arg(numberOfControlPointsPerRow).arg(order).arg(derivativeU).arg(derivativeV).arg(u0).arg(u1).arg(v0).arg(v1).arg(paramA).arg(paramB);
    //const QString filename = m_defaultFolder + QString("/energymatrix.mat.dm");

    const int order = properties->getOrder();
    const int numberOfControlPointsPerRow = properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerPatch = numberOfControlPointsPerRow * numberOfControlPointsPerRow;

    const QString filename = this->createFilenameSurfaceUniformKnots(properties, derivativeU, derivativeV, u0, u1, v0, v1);
    if (m_loadFromFile && QFile(filename).exists())
        return this->loadSparseMatrixFromFile(filename, numberOfControlPointsPerPatch, numberOfControlPointsPerPatch);



    const QVector<double> *knots = properties->getKnotsP();
    const int lastInnerKnotIndex = properties->getLastInnerKnotIndex();

    const int numberOfGaussPointsU = order - derivativeU + 5;   //only +1 needed, but we want to prevent numerical errors
    const int numberOfGaussPointsV = order - derivativeV + 5;
    gsl_integration_glfixed_table *tableU = gsl_integration_glfixed_table_alloc(numberOfGaussPointsU);
    gsl_integration_glfixed_table *tableV = gsl_integration_glfixed_table_alloc(numberOfGaussPointsV);

    Eigen::SparseMatrix<double> matrix(numberOfControlPointsPerPatch, numberOfControlPointsPerPatch);
    Eigen::MatrixXd matrixDense = Eigen::MatrixXd::Zero(numberOfControlPointsPerPatch, numberOfControlPointsPerPatch);
    Eigen::VectorXi nonZeros = Eigen::VectorXi::Zero(numberOfControlPointsPerPatch);

    for (int i1j1 = 0; i1j1 < numberOfControlPointsPerPatch; i1j1++) {
        int i1, j1;
        BSplineSurface::vecToMatIndexLocal(i1j1, i1, j1, numberOfControlPointsPerRow);
        for (int i2j2 = i1j1; i2j2 < numberOfControlPointsPerPatch; i2j2++) {
            int i2, j2;
            BSplineSurface::vecToMatIndexLocal(i2j2, i2, j2, numberOfControlPointsPerRow);

            double integralU = 0;
            double integralV = 0;
            //for (int iKnots = 0; iKnots < sizeOfKnotVector - 1; iKnots++) {
            for (int iKnots = properties->getFirstInnerKnotIndex(); iKnots < lastInnerKnotIndex; iKnots++) { //only consider inner knot vectors
                const double a = knots->at(iKnots);
                const double b = knots->at(iKnots + 1);
                if (a != b) {
                    double ua = a;
                    double ub = b;
                    if (ua < u0)
                        ua = u0;
                    if (ub > u1)
                        ub = u1;
                    if (ub > ua) {
                        for (int i = 0; i < numberOfGaussPointsU; i++) {
                            double xi, wi;
                            gsl_integration_glfixed_point(ua, ub, i, &xi, &wi, tableU);
                            integralU += wi * properties->NDerivative(xi, i1, derivativeU) * properties->NDerivative(xi, i2, derivativeU);
                        }
                    }

                    double va = a;
                    double vb = b;
                    if (va < v0)
                        va = v0;
                    if (vb > v1)
                        vb = v1;
                    if (vb > va) {
                        for (int i = 0; i < numberOfGaussPointsV; i++) {
                            double xi, wi;
                            gsl_integration_glfixed_point(va, vb, i, &xi, &wi, tableV);
                            integralV += wi * properties->NDerivative(xi, j1, derivativeV) * properties->NDerivative(xi, j2, derivativeV);
                        }
                    }
                }
            }

            if (integralU * integralV != 0) {
                //matrix.insert(i1j1, i2j2) = integralU * integralV;
                matrixDense(i1j1, i2j2) = integralU * integralV;
                nonZeros(i1j1) += 1;
                if (i1j1 != i2j2) {
                    matrixDense(i2j2, i1j1) = integralU * integralV;
                    nonZeros(i2j2) += 1;
                    //matrix.insert(i2j2, i1j1) = integralU * integralV;
                }
            }
        }
    }

    matrix.reserve(nonZeros);
    for (int i = 0; i < numberOfControlPointsPerPatch; i++)
        for (int j = 0; j < numberOfControlPointsPerPatch; j++) {
            double val = matrixDense(i, j);
            if (val != 0)
                matrix.insert(i, j) = matrixDense(i, j);
        }

    gsl_integration_glfixed_table_free(tableU);
    gsl_integration_glfixed_table_free(tableV);

    //TODO check if dir exists first (and create it if it does not exist)
    if (m_saveToFile) {
        if (!QFile(filename).exists())
            this->saveSparseMatrixToFile(filename, matrix);
    }

//    double offset = 0.0000000000001;
//    qDebug() << "Added offset to energy matrix" << offset;
//    for (int i = 0; i < numberOfControlPointsPerPatch; i++)
//        matrix.coeffRef(i, i) += offset;

//    bool eigenValuesOk = false;
//    while (!eigenValuesOk) {
//        eigenValuesOk = true;
//        Eigen::EigenSolver<Eigen::MatrixXd> eSolver;
//        eSolver.compute(m_energyMatrix.toDense(), true);
//        Eigen::VectorXcd eigenvalues = eSolver.eigenvalues();
//        for (int i = 0; i < m_numberOfControlPointsPerPatch; i++)
//            if (eigenvalues[i].real() < 0) {
//                Eigen::VectorXcd eigenVec = eSolver.eigenvectors().col(i);
//                //energy + eigenVec^T * v has eigenValue eVal + v^T * eigenvec
//            }
//        std::cout << "Eigenvalues" << std::endl << eigenvalues << std::endl;
//    }
//    qDebug() << "rank one update on energy matrix needs to be finished and tested";

    return matrix;
}

QString EnergyMatrixGenerator::createFilenameSurfaceUniformKnots(BSplineProperties *properties, const int derivativeU, const int derivativeV, const double u0, const double u1, const double v0, const double v1)
{
    const int order = properties->getOrder();
    const int nc = properties->getNumberOfControlPoints();
    const double knotStart = properties->getInnerKnotsStart();
    const double knotEnd = properties->getInnerKnotsEnd();

    return QString(m_defaultFolder + "/surface/uniform_s%1_e%2/u%3_%4_v%5_%6/matrix_o%7_nc%8_ndu%9_ndv_%10.bin").arg(knotStart).arg(knotEnd).arg(u0).arg(u1).arg(v0).arg(v1).arg(order).arg(nc).arg(derivativeU).arg(derivativeV);
}

void EnergyMatrixGenerator::saveSparseMatrixToFile(QString filename, Eigen::SparseMatrix<double> &matrix)
{
    const QString dirPath(QDir::cleanPath(filename + "/.."));
    if (!QDir(dirPath).exists())
        QDir().mkpath(dirPath);

    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);

    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            const double value = matrix.coeff(i, j);
            QByteArray array(reinterpret_cast<const char*>(&value), sizeof(double));
            file.write(array);
        }
    }

    file.close();
}

Eigen::SparseMatrix<double> EnergyMatrixGenerator::loadSparseMatrixFromFile(QString filename, const int numberOfRows, const int numberOfCols)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly))
        return Eigen::SparseMatrix<double>(1, 1);

    Eigen::SparseMatrix<double> matrix(numberOfRows, numberOfCols);
    const int doubleSize = sizeof(double);
    const QByteArray byteNum = file.read(numberOfRows * numberOfCols * doubleSize);

    int position = 0;
    for (int i = 0; i < numberOfRows; i++) {
        for (int j = 0; j < numberOfCols; j++) {
            const double value = *reinterpret_cast<const double*>(byteNum.data() + position);
            position += doubleSize;

            if (value != 0)
                matrix.insert(i, j) = value;
        }
    }

    file.close();

    return matrix;
}

#ifndef ENERGYMATRIXGENERATOR_H
#define ENERGYMATRIXGENERATOR_H

#include "BSpline/bsplineproperties.h"

//Eigen includes
#include <SparseCore>

//Qt incluces
#include <QPair>
#include <QString>

class EnergyMatrixGenerator
{
public:
    EnergyMatrixGenerator();
    EnergyMatrixGenerator(const QString &defaultFolder, const bool loadFromFileIfPossible = true, const bool saveToFileIfPossible = true);
    //TODO add setter for load/safe and default folder

    //TODO add load/safe option for curve case
    Eigen::SparseMatrix<double> getEnergyMatrixCurve(BSplineProperties *properties, const int derivative, const double t0, const double t1);
    Eigen::SparseMatrix<double> getEnergyMatrixSurface(BSplineProperties *properties, const int derivativeU, const int derivativeV, const double u0, const double u1, const double v0, const double v1);
private:
    QString createFilenameSurfaceUniformKnots(BSplineProperties *properties, const int derivativeU, const int derivativeV, const double u0, const double u1, const double v0, const double v1);
    void saveSparseMatrixToFile(QString filename, Eigen::SparseMatrix<double> &matrix);
    Eigen::SparseMatrix<double> loadSparseMatrixFromFile(QString filename, const int numberOfRows, const int numberOfCols);

    QString m_defaultFolder;
    bool m_loadFromFile;
    bool m_saveToFile;
};

#endif // ENERGYMATRIXGENERATOR_H

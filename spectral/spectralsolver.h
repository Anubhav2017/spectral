#ifndef SPECTRALSOLVER_H
#define SPECTRALSOLVER_H
#include<laplacesolver.h>
#include<fstream>

class spectralsolver
{
public:

    static QVector<double> computeSpectralEigenvector(Mesh *m, int ef);
    static void saveEigenvectorToTxt(Mesh *m, QVector<double> eigenvector, QString filename);
    static QVector<int> loadMorseSmaleDecomposition(QString filename);
    static CellMesh* createCellMeshFromMorseSmale(Mesh *m, QVector<int> vertexmap);
    static QVector<QColor> generateColorMap2(QVector<double >realvector);
};

#endif // SPECTRALSOLVER_H

#ifndef SPLINEEVALUATOR_H
#define SPLINEEVALUATOR_H

#include "BSpline/bsplinepatchnetwork.h"
#include "CellMesh/cellmesh.h"
#include "Mesh/mesh.h"

#include <QVector>
#include <SparseCore>

class SplineEvaluator
{
public:
    SplineEvaluator(BSplinePatchNetwork *bSplinePatchNetwork, Eigen::SparseMatrix<double> *energyMatrix, int numberOfEvalPointsForG1, bool doOptimizationForDataError);

    void evaluateLSErrors(bool doOptimizationForDataError);
    void evaluateG1Errors(int numberOfG1SamplePoints);
    void evaluatePatchEnergyValues();
    void evaluatePatchAreas();

    void writeQualityFile(QString filename);
    void writeMatlabTriangulationStructureFile(QString filename);
    void writeMatlabRegularSamplingStructureFile(int samplePointsPerRow, QString filename);
    void writeMatlabRegularSampledCurvature(int samplePointsPerRow, QString filename);
    void writeMatlabRegularSamplingG1Error(QString filename);
    void writeMatlabDataErrorFiles(QString filenameTri, QString filenameBSpline);
    void writeMatlabG1ErrorAtPLVFiles(QString filenameTemplateTri, QString filenameTemplateBSpline);

    double getRmsError();
    double getMaxDataError();
    double getTotalAverageDetError();
    double getTotalAverageNormalError();
    double getTotalAverageAngleError();
    double getMaxDetError();
    double getMaxNormalError();
    double getMaxAngleError();
    double getBoundingBoxScaledRmsError();
    double getBoundingBoxScaledMaxDataError();
    double getBoundingSphereScaledRmsError();
    double getBoundingSphereScaledMaxDataError();
    double getTotalSurfaceAreaMesh();
    double getTotalSurfaceAreaBSplines();
    double getTotalEnergy();

    static double calculatePatchEnergy(BSplineSurface *surface, Eigen::SparseMatrix<double> &energyMatrix);
private:
    BSplinePatchNetwork *m_bSplinePatchNetwork;
    CellMesh *m_cellMesh;
    Mesh *m_mesh;
    Eigen::SparseMatrix<double> *m_energyMatrix;

    double m_meshSizeX;
    double m_meshSizeY;
    double m_meshSizeZ;

    QVector3D m_boundingSphereCenter;
    double m_boundingSphereRadius;

    QVector<double> m_dataErrors;
    QVector<QVector3D> m_projectedCoordinates;
    double m_totalSquaredError;
    double m_rmsError;
    double m_maxDataError;
    double m_averageEdgeLength;
    double m_boundingBoxScaledRmsError;
    double m_boundingBoxScaledMaxDataError;
    double m_boundingSphereScaledRmsErrorPercent;
    double m_boundingSphereScaledMaxDataErrorPercent;

    QVector<double> m_patchEnergyValues;
    QVector<double> m_patchAreasMesh;
    QVector<double> m_patchAreasBSplines;

    double m_totalEnergy;
    double m_totalAreaMesh;
    double m_totalAreaBSplines;

    QVector<QVector<double> > m_individualDetErrors;
    QVector<QVector<double> > m_individualNormalErrors;
    QVector<QVector<double> > m_individualAngleErrors;
    QVector<QVector<QVector3D> > m_individualG1ErrorPositions;
    double m_maxDetError;
    double m_maxNormalError;
    double m_maxAngleError;
    QVector<double> m_integratedDetErrors;
    QVector<double> m_integratedNormalErrors;
    QVector<double> m_integratedAngleErrors;
    QVector<double> m_individualBorderCurveLengths;
    double m_totalIntegratedDetError;
    double m_totalIntegratedNormalError;
    double m_totalIntegratedAngleError;
    double m_totalBorderCurveLength;
    double m_totalAverageDetError;
    double m_totalAverageNormalError;
    double m_totalAverageAngleError;

};

#endif // SPLINEEVALUATOR_H

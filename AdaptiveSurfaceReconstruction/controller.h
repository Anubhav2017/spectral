#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "Mesh/mesh.h"
#include "CellMesh/cellmesh.h"
#include "CellMesh/parameterization.h"
#include "BSpline/bsplinesurface.h"
#include "BSpline/bsplinecurve.h"
#include "BSpline/bsplinepatchnetwork.h"
#include "BSpline/bsplineproperties.h"
#include "MBSI/subdivisionabstractmetric.h"

#include "stlwriter.h"

#include <SparseCore>

#include <QVector>
#include <QString>
#include <QTime>
#include <QVector3D>
#include <QPair>
#include <QColor>

class Controller
{
public:
    enum QuadmeshAlgorithm {QMA_DelaunayMatching   = 0,
                            QMA_VoronoiSubdivision = 1,
                            QMA_CuttingPlanes      = 2};

    enum BSplineAlgorithm {BSA_CurveFittingG0          = 0,  //cg0
                           BSA_CurveFittingG0Iterative = 1,  //cg0i
                           BSA_GlobalG0                = 2,  //gg0
                           BSA_SimplifiedGlobalG1      = 3,  //gsg1
                           BSA_LocalG1Multiphase       = 4}; //lg1

    Controller();
    ~Controller();

    void start();

    void setStlWriter(StlWriter *writer);
    void setDebugingMode(int mode);
    void setConsoleOutput(bool mode);
    void setLogOutput(bool mode);
    void setStartPoint(int startPoint);
    void setEndPoint(int endPoint);
    void setDmFolder(QString dmFolder);
    void setEnergyFolder(QString energyFolder);
    void setStlBSplineSamplePoints(int samplePoints);
    void setStlLineWidth(double lineWidth);
    void setNumberOfThreads(int numberOfThreads);

    void setParameterInterval(double a, double b);

    void setMeshScaling(int scaling);
    void setKeepingProportions(bool keepProportions);
    void setScalingFactor(double scalingFactor);

    void setmaxSmoothIterations(int maxIterations);

    void setQuadmeshAlgorithm(QuadmeshAlgorithm alg);
    void setMinNumberOfInitialCells(int numberOfCells);
    void setMinNumberOfQuadCells(int numberOfCells);
    void setInstantTopologyCheck(bool instantTopologyCheck);
    void setSkipRefinement(bool skipRefinement);
    void setNumberOfCentralSubdivisions(int numberOfSubdivisions);

    void setNumberOfVertexOptimizations(int numberOfOptimizations);

    void setBSplineAlgorithm(BSplineAlgorithm alg);
    void setRecomputeParameterization(bool recomputeParam);
    void setInvertCellmeshOrientation(bool invert);
    void setUsageOfInitialGuess(bool useInitialGuess);
    void setNumberOfControlPoints(int numberOfControlPointsPerRow);
    void setSplineOrder(int splineOrder);
    void setWeightEnergyTerm(double weight);
    void setG1VertexParameterRange(double range);
    void setNumberOfG1VertexPoints(int numberOfPoints);
    void setVertexFittingG1Opt(bool mode);
    void setVertexFitG1Weight(double weight);
    void setNumberOfG1EdgePoints(int numberOfPoints);
    void setEdgeFittingG1Opt(bool mode);
    void setEdgeFitG1Weight(double weight);
    void setDetailedFittingOutput(bool detailedFittingOutput);
    void setEnergyMatrixStorageOnHarddrive(bool storeMatrices);
    void setEnergyMatrixLoadingFromHarddrive(bool loadMatrices);

    void setNumberOfParameterOptimizations(int numberOfOptimizations);
    void setUseOfNLOforParamOpt(bool useNLO);

    void setScalingForEvaluation(int scaling);
    void setMatlabOutput(bool doMatlabOutput);
    void set3DImageOutput(bool do3DImage);
    void set3DImageSizes(int sizeX, int sizeY, int sizeZ);

    void setMbsiMetricType(int metricType);
    void setMbsiThreshold(double thresholdFactor, bool factorOfAverage);
    void setMbsiTargetNumberOfVPC(int number);
    void setMbsiDoReferenceAlgorithms(bool doAlgorithms);
    void setMbsiPivotPointSize(double size);
    void setMbsiCamPointSize(double size);
    void setMbsiCamDirectionSize(double size);
    void setMbsiCamDirectionLength(double length);
    void setMbsiCamDistance(double distance);

    void setNloptMaxIterations(int maxIterations);
    void setNloptFPrec(double fPrec);
    void setNloptXPrec(double xPrec);
    void setNloptConstrPrec(double constrPrec);

    static bool loadBasepointsAndNormalsFromFile(QString filename, QVector<QVector3D> *basepoints, QVector<QVector3D> *normals);

    void createMatlab2DOutputForMetric(BSplineSurface *surf, QString csvFilename, SubdivisionAbstractMetric *metric, SubSurfaceTree *subdivision, const int numberOfEvalPoints, QPair<double, double> parameterInterval);

    void writeOutput(const QString msg) const;
    void writeToConsoleOnly(const QString msg) const;
    void writeToLogOnly(const QString msg) const;
    void clearLogFile();

    static bool removeRecursively(const QString &dirName);
private:
//    void createDebugSpline();

    void cellMeshToStl(QString basename, QTime *timeStep, bool storeIndividualCellVerticesAndEdges, bool storeCellContents);
    void bSplinesToStl(int numberOfEvalPoints, QString basename, QTime *timeStep, bool doControlPoints, bool doKnotLines, bool doDataError, bool doDataErrorForIndividualSurfaces);

    StlWriter *m_stlWriter;
    bool m_consoleOutput;
    bool m_logOutput;
    QString m_logFilename;
    Mesh *m_mesh;
    CellMesh *m_cellMesh;
    BSplinePatchNetwork *m_bSplinePatchNetwork;

    int m_startPoint;
    int m_endPoint;
    int m_debugingMode;
    QString m_dmFolder;
    QString m_energyFolder;
    int m_stlBSplineSamplePoints;
    double m_stlLineWidth;
    int m_fixedNumberOfThreads;

    bool m_useLoadedParamInterval;
    QPair<double, double> m_parameterInterval;

    int m_scaleMesh;
    bool m_keepMeshProportions;
    double m_scalingFactor;

    int m_maxSmoothIterations;

    QuadmeshAlgorithm m_quadmeshAlgorithm;
    int m_minNumberOfInitialCells;
    int m_minNumberOfQuadCells;
    int m_numberOfCentralSubdivisions;
    bool m_instantTopologyCheck;
    bool m_skipRefinement;

    int m_numberOfVertexOptimizations;

    BSplineAlgorithm m_bSplineAlgorithm;
    BSplineProperties *m_bSplineProperties;
    bool m_recomputeParameterization;
    bool m_invertCellmeshOrientation;
    bool m_useInitialGuess;
    int m_numberOfControlPoints;
    int m_splineOrder;
    double m_weightEnergyTerm;
    int m_numberOfG1EdgePoints;
    bool m_g1VertexParameterRangeSet;
    double m_g1VertexParameterRange;
    int m_numberOfG1VertexPoints;
    bool m_doVertexFittingWithG1Opt;
    double m_vertexFitG1Weight;
    bool m_doEdgeFittingWithG1Opt;
    double m_edgeFitG1Weight;
    bool m_detailedFittingOutput;
    bool m_storeEnergyMatricesOnHardDrive;
    bool m_loadEnergyMatricesIfPossible;

    int m_numberOfParameterOptimizations;
    bool m_useNLOforParameterOpt;

    int m_scaleForEvaluation;
    bool m_doMatlabOutput;
    bool m_do3DImageOutput;
    int m_3DimageSizeX;
    int m_3DimageSizeY;
    int m_3DimageSizeZ;

    int m_mbsiMetricType;
    double m_mbsiThreshold;
    int m_mbsiTargetNumberOfVPC;
    bool m_mbsiThresholdIsAVGFactor;
    bool m_mbsiDoReferenceAlgorithms;
    double m_mbsiPivotPointSize;
    double m_mbsiCamPointSize;
    double m_mbsiCamDirVectorSize;
    double m_mbsiCamDirVectorLength;
    double m_mbsiCamDistance;

    int m_nloptMaxIterations;
    double m_nloptFPrec;
    double m_nloptXPrec;
    double m_nloptConstrPrec;
};

#endif // CONTROLLER_H

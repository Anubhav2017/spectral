#include "controller.h"

#include "SurfaceTessellation/surfacevoronoigenerator.h"
#include "SurfaceTessellation/quadmeshgenerator.h"
#include "BSpline/bsplinetriangulation.h"
#include "MBSI/viewpointcandidategenerator.h"
#include "MBSI/subdivisionthinplateenergymetric.h"
#include "MBSI/subdivisionexactcurvaturemetric.h"
#include "MBSI/subdivisionsurfaceareametric.h"
#include "MBSI/subdivisionmanualhighlightmetric.h"
#include "MBSI/subdivisionnormaldeviationmetric.h"
#include "MBSI/subdivisioninspectionanglemetric.h"
#include "MBSI/colormapgenerator.h"
#include "MBSI/surfacesubdivider.h"
#include "parameteroptimizer.h"
#include "energymatrixgenerator.h"
#include "splinefitter.h"
#include "splineevaluator.h"
#include "voxelizer.h"
#include "plywriter.h"

#include <iostream>
#include <limits>

#include <QQueue>

#include <QDir>
#include <QFileDialog>
#include <QDebug>

#include "IterativeLinearSolvers"   //includes least squares solver which is only dev. version right now
#include "SparseCholesky"
#include "SparseLU"
#include <Sparse>
#include <Dense>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <omp.h>

#include <nlopt.h>

#include "external/Miniball.hpp"
#include "external/DataIO.h"

Controller::Controller():
    m_stlWriter(0), m_mesh(0), m_cellMesh(0), m_bSplinePatchNetwork(0), m_bSplineProperties(0)
{
    // individual steps
    // 0 - load isosurface from stl or compute it from data
    // 1 - do mesh smoothing
    // 2 - compute the quadmesh
    // 3 - optimize the quadmesh
    // 4 - do BSpline Approximation
    // 5 - optimize parameterization by projection onto tangential plane
    // 6 - compute quality of the b-splines
    // 7 - model based surface inspection stuff


    m_startPoint = 0;
    m_endPoint = 6;
    m_debugingMode = 1;
    m_dmFolder = QString("../dm/");
    m_energyFolder = QString("../energy/");
    m_stlBSplineSamplePoints = 51;
    m_stlLineWidth = 0.3;
    m_fixedNumberOfThreads = -1;    //set to < 0 to let the system choose

    m_consoleOutput = true;
    m_logOutput = true;
    m_logFilename = QString(m_dmFolder + "log.txt.dm"); //this also gets changed when the dmFolder gets changed

    m_useLoadedParamInterval = true;
    m_parameterInterval = QPair<double, double>(0, 1);

    //0 - Isosurface parameters
    m_scaleMesh = -1;
    m_keepMeshProportions = false;
    m_scalingFactor = 1;

    //1 - Smoothing parameters
//    const bool doMeshSmoothing = false;
    m_maxSmoothIterations = 0;
//    const double maxSmoothDistFromOrigin = 1;
//    const double abortSmoothAverageDist = 0.01;

    //2 - Initial quadmesh parameters
//    Parameterization::WeightType param_wt = Parameterization::WT_cotangent;
//    Parameterization::BorderType param_bt = Parameterization::BT_lengthWeightedUnitSquare;
    m_quadmeshAlgorithm = QMA_DelaunayMatching;
    m_minNumberOfInitialCells = 6;
    m_minNumberOfQuadCells = 6;
    m_instantTopologyCheck = false;
    m_numberOfCentralSubdivisions = 0;
    m_skipRefinement = false;

    //3 - Quadmesh optimization parameters
//    const bool doQuadmeshOptimization = true;
    m_numberOfVertexOptimizations = 1;

    //4 - B-Spline parameters
    m_bSplineAlgorithm = BSA_LocalG1Multiphase;
    m_recomputeParameterization = false;
    m_invertCellmeshOrientation = false;
    m_useInitialGuess = false;
    m_numberOfControlPoints = 6;
    m_splineOrder = 3;
    m_weightEnergyTerm = 0.5;
    //vertex fit parameters
    m_g1VertexParameterRangeSet = false;
    m_g1VertexParameterRange = 0.25;
    m_numberOfG1VertexPoints = 4;
    m_doVertexFittingWithG1Opt = false;
    m_vertexFitG1Weight = 1000;
    //edge fit parameters
    m_numberOfG1EdgePoints = m_numberOfControlPoints - 4;
    m_doEdgeFittingWithG1Opt = false;
    m_edgeFitG1Weight = 1000;
    //other
    m_detailedFittingOutput = false;
    //energy matrix
    m_storeEnergyMatricesOnHardDrive = false;
    m_loadEnergyMatricesIfPossible = true;

    //5 - Optimize parameters (projection onto tangential space)
    m_useNLOforParameterOpt = false;
    m_numberOfParameterOptimizations = 1;

    //6 - Quality measurement parameters
//    const int numberOfEvalPointsForG1 = 101;
//    const int numberOfEvalPointsForCurvature = 51;
//    const bool doOptimizationForDataError = true;
    m_scaleForEvaluation = -1;
    m_doMatlabOutput = true;
    m_do3DImageOutput = false;
    m_3DimageSizeX = -1;
    m_3DimageSizeY = -1;
    m_3DimageSizeZ = -1;

    //7 - Model based surface inspection
    m_mbsiMetricType = 0;
    m_mbsiThreshold = 0.5; //subdivision threshold = factor * energy / numberOfPatches
    m_mbsiThresholdIsAVGFactor = true;
    m_mbsiTargetNumberOfVPC = -1;
    m_mbsiDoReferenceAlgorithms = false;
    m_mbsiPivotPointSize = 0.6;
    m_mbsiCamPointSize = 0.4;
    m_mbsiCamDirVectorSize = 0.2;
    m_mbsiCamDirVectorLength = 0;
    m_mbsiCamDistance = 2;

    //Parameters for numerics
    m_nloptMaxIterations = 1000000;
    m_nloptFPrec = 0.0000001;
    m_nloptXPrec = 0.0000001;
    m_nloptConstrPrec = 0.001;

}

Controller::~Controller()
{
    if (m_bSplinePatchNetwork)
        delete m_bSplinePatchNetwork;
    if (m_bSplineProperties)
        delete m_bSplineProperties;
    if (m_cellMesh)
        delete m_cellMesh;
    if (m_mesh)
        delete m_mesh;
}


void Controller::start()
{
    //Start timer
    QTime timeTotal, timeStep;
    timeTotal.start();
    timeStep.start();

    //Hardcoded variables (can be overwriten by user parameters)
    const double maxPercentageOfThreads = 0.8;

    //0 - Isosurface parameters

    //1 - Smoothing parameters
    const double maxSmoothDistFromOrigin = 5;
    const double abortSmoothAverageDist = 0.001;
    const bool doMeshSmoothing = (m_maxSmoothIterations > 0);

    //2 - Initial quadmesh parameters
    Parameterization::WeightType param_wt = Parameterization::WT_cotangent;
    Parameterization::BorderType param_bt_forQuadmeshGeneration = Parameterization::BT_lengthWeightedCircle;
    Parameterization::BorderType param_bt_forSplineFitting = Parameterization::BT_lengthWeightedUnitSquare;
    int parameterizationStorageReservation_ExpectedEdgesPerVertex = -1; //number of expected edges per vertex in the mesh (gets initialized with smart default values in cell mesh as well)

    //3 - Quadmesh optimization parameters

    //4 - B-Spline parameters
    m_bSplineProperties = BSplineProperties::createInstanceUniformKnotVectorMultipleEnds(m_splineOrder, m_numberOfControlPoints,
                                                                                         m_parameterInterval.first, m_parameterInterval.second);
    QVector<QPair<int, int> > energyDerivatives;
    QVector<int> energyDerivativeFactors;
    energyDerivatives << QPair<int, int>(2, 0) << QPair<int, int>(1, 1) << QPair<int, int>(0, 2);
    energyDerivativeFactors << 1 << 2 << 1;
//    energyDerivatives << QPair<int, int>(3, 0) << QPair<int, int>(0, 3);
//    energyDerivativeFactors << 1 << 1;

    //5 - Optimize parameters (projection onto tangential space)
    const bool doParamOpt = true;

    //6 - Quality measurement parameters
    const int numberOfEvalPointsForG1 = 101;
    const int numberOfEvalPointsForCurvature = 101;
    const bool doOptimizationForDataError = true;
    const int iassSurfaceSamplePoints = 101;
    const int iassPadding = 5;

    //7 - Model based surface inspection
    if (m_mbsiCamDirVectorLength == 0)
        m_mbsiCamDirVectorLength = m_mbsiCamDistance;

    //DEBUG parameters
    //If no stlWriter is set, deactivate debuging mode
    if (!m_stlWriter) {
        qDebug() << "No stl writer set. Debuging mode set to 0!";
        m_debugingMode = 0;
    }

    //Set individual outputs according to debuging mode
    const bool storeImportedMeshAsStl = (m_debugingMode > 0);
    const bool storeSmoothedMeshAsStl = (m_debugingMode > 0);

    const bool storeImportedCellMeshAsStl = (m_debugingMode > 0);
    const bool storeInitialCellMeshAsStl = (m_debugingMode > 0);
    const bool storeOptimizedCellMeshAsStl = (m_debugingMode > 0);
    const bool cellMeshIndividualVerticesAndEdgesToStl = (m_debugingMode > 1);
    const bool cellMeshContentsToStl = (m_debugingMode > 2);

    const bool saveParameterizationsAsImage = (m_debugingMode > 1);

    const bool storeImportedBSplinesAsStl = (m_debugingMode > 0);
    const bool storeComputedBSplinesAsStl = (m_debugingMode > 0);
    const bool bSplineControlPointsToStl = (m_debugingMode > 0);
    const bool bSplineKnotLinesToStl = (m_debugingMode > 1);
    const bool bSplineDataErrorToStl = (m_debugingMode > 1);
    const bool bSplineDataErrorOfIndividualSurfacesToStl = (m_debugingMode > 2);

    const bool writeMBSISubdivisionsToSTL = (m_debugingMode > 0);
    const bool writeMBSIViewpointCandidatesToSTL = (m_debugingMode > 0);
    const bool doExtendedMBSIOutput = (m_debugingMode > 1);

    //Other parameters
    const bool exportMesh = true;
    const bool exportCellMesh = true;
    const bool exportParameterization = true;
    const bool exportBSplines = true;

    const bool checkCellMeshAfterLoading = true;

    //Filenames for import and export
    const QString paramSubfolder("parameterizations/");
    const QString matlabSubfolder("matlab/");
    const QString meshOpenFilename(m_dmFolder + "mesh.m.dm");
    const QString meshSaveFilename(m_dmFolder + "mesh.m.dm");
    const QString cellMeshOpenFilename(m_dmFolder + "cellmesh.cm.dm");
    const QString cellMeshSaveFilename(m_dmFolder + "cellmesh.cm.dm");
    const QString cellMeshBackupSaveFilename(m_dmFolder + "cellmesh_backup.cm.dm");
    const QString bSplineOpenFilename(m_dmFolder + "bsplines.bs.dm");
    const QString bSplineSaveFilename(m_dmFolder + "bsplines.bs.dm");
    const QString bSplineBackupSaveFilename(m_dmFolder + "bsplines_backup.bs.dm");
    const QString parameterizationOpenFilename(m_dmFolder + paramSubfolder + "param_%1.p.dm");
    const QString parameterizationSaveFilename(m_dmFolder + paramSubfolder + "param_%1.p.dm");
    const QString parameterizationImageFilename(m_dmFolder + paramSubfolder + "paramExport_%1.png");
    const QString qualityFilename(m_dmFolder + "quality.txt.dm");
    const QString matlabTriStructureFilename(m_dmFolder + matlabSubfolder + "triangulation.csv");
    const QString matlabTriErrorFilename(m_dmFolder + matlabSubfolder + "triError.csv");
    const QString matlabBSErrorFilename(m_dmFolder + matlabSubfolder + "bsplineError.csv");
    const QString matlabTriBordersFilename(m_dmFolder + matlabSubfolder + "triBorders_%1.csv");
    const QString matlabBSBorderFilename(m_dmFolder + matlabSubfolder + "bsplineBorders_%1.csv");
    const QString matlabRegularSamplingBorderFilename(m_dmFolder + matlabSubfolder + "regularSamplingBorders.csv");
    const QString matlabRegularSamplingStructureFilename(m_dmFolder + matlabSubfolder + "regularSamplingStructure.csv");
    const QString matlabRegularSamplingCurvatureFilename(m_dmFolder + matlabSubfolder + "regularSamplingCurvature.csv");
    const QString iassOutputFilename(m_dmFolder + "result3DImage.iass");

    if (!QDir(m_energyFolder).exists())
        QDir().mkpath(m_energyFolder);

    if (!QDir(m_dmFolder).exists())
        QDir().mkpath(m_dmFolder);

    if (m_recomputeParameterization)
        Controller::removeRecursively(m_dmFolder + paramSubfolder);

    if (!QDir(m_dmFolder + paramSubfolder).exists())
        QDir().mkpath(m_dmFolder + paramSubfolder);

    if (!QDir(m_dmFolder + matlabSubfolder).exists())
        QDir().mkpath(m_dmFolder + matlabSubfolder);

    //Setup number of threads for omp
    const int numberOfCores = omp_get_num_procs();
    int numberOfThreads = floor(numberOfCores * maxPercentageOfThreads);
    if (numberOfThreads < 1)
        numberOfThreads = 1;
    if (m_fixedNumberOfThreads > 0)
        numberOfThreads = m_fixedNumberOfThreads;
    omp_set_num_threads(numberOfThreads);

    writeOutput(QString("Process started! Number of threads: %1").arg(numberOfThreads));

    //Check if parameters make sense
    const double recommendedRange = m_bSplineProperties->getKnotsP()->at(m_splineOrder + 2);
    if (!m_g1VertexParameterRangeSet && m_startPoint <= 4 && m_endPoint >= 4 && m_g1VertexParameterRange < recommendedRange)
        m_g1VertexParameterRange = recommendedRange;

    //import necessary data
    if (m_startPoint > 0) {
        m_mesh = Mesh::loadFromFile(meshOpenFilename);
        if (!m_mesh) {
            writeOutput(QString("ERROR: Could not open mesh file: %1").arg(meshOpenFilename));
            exit(-1);
        }
        writeOutput(QString("mesh imported. Time: %1 ms").arg(timeStep.restart()));
        writeOutput(QString("Vertices: %1, Edges: %2, Faces: %3").arg(m_mesh->getNumberOfVertices()).arg(m_mesh->getNumberOfEdges()).arg(m_mesh->getNumberOfFaces()));

        if (storeImportedMeshAsStl) {
            m_stlWriter->writeMeshFacesToStl(m_mesh, "importedMesh.stl");
            writeOutput(QString("imported Mesh stored. Time: %1 ms").arg(timeStep.restart()));
        }
    }
    if (m_startPoint > 2) {
        m_cellMesh = CellMesh::loadFromFile(cellMeshOpenFilename, m_mesh);
        if (!m_cellMesh) {
            writeOutput(QString("ERROR: Could not open cellmesh file: %1").arg(cellMeshOpenFilename));
            exit(-1);
        }
        if (m_useLoadedParamInterval)
            m_parameterInterval = m_cellMesh->getParameterInterval();
        else
            m_cellMesh->setParameterInterval(m_parameterInterval);

        if (parameterizationStorageReservation_ExpectedEdgesPerVertex > 0)
            m_cellMesh->setParameterizationStorageReservation(parameterizationStorageReservation_ExpectedEdgesPerVertex);
        writeOutput(QString("cellMesh imported. Time: %1 ms").arg(timeStep.restart()));

        if (checkCellMeshAfterLoading) {
            m_cellMesh->checkMesh();
            writeOutput(QString("cellMesh checked. Time: %1 ms").arg(timeStep.restart()));
        }

        if (!m_recomputeParameterization) {
            m_cellMesh->loadParameterizations(parameterizationOpenFilename);
            writeOutput(QString("parameterizations imported. Time: %1 ms").arg(timeStep.restart()));
        } else {
            m_cellMesh->computeFullParameterization(param_bt_forQuadmeshGeneration);
        }

        if (storeImportedCellMeshAsStl) {
            this->cellMeshToStl(QString("importedCellMesh"), &timeStep, cellMeshIndividualVerticesAndEdgesToStl, cellMeshContentsToStl);
            writeOutput(QString("imported cellmesh stored. Time: %1 ms").arg(timeStep.restart()));
        }
    }
    if (m_startPoint > 4 || (m_useInitialGuess && m_endPoint >= 4)) {
        m_bSplinePatchNetwork = BSplinePatchNetwork::loadFromFile(bSplineOpenFilename, m_cellMesh);
        if (!m_bSplinePatchNetwork) {
            writeOutput(QString("ERROR: Could not open B-spline file: %1").arg(bSplineOpenFilename));
            exit(-1);
        }

        if (m_bSplineProperties)
            delete m_bSplineProperties;

        m_bSplineProperties = m_bSplinePatchNetwork->getProperties();
        m_numberOfControlPoints = m_bSplineProperties->getNumberOfControlPoints();
        m_splineOrder = m_bSplineProperties->getOrder();

        writeOutput(QString("bSplines loaded. Number of CP per row: %1, order: %2. Time: %3 ms").arg(m_numberOfControlPoints).arg(m_splineOrder).arg(timeStep.restart()));

        if (storeImportedBSplinesAsStl)
            this->bSplinesToStl(m_stlBSplineSamplePoints, QString("importedBSpline"), &timeStep, bSplineControlPointsToStl, bSplineKnotLinesToStl, bSplineDataErrorToStl, bSplineDataErrorOfIndividualSurfacesToStl);
    }

    //0 - Load or compute the isosurface
    if (m_startPoint <= 0 && m_endPoint >= 0) {
        QString filename = QFileDialog::getOpenFileName(0, QString("Load triangulation"), QString("../.."), "*.stl");
        if (filename.isEmpty()) {
            writeOutput(QString("No .stl file selected"));
            exit(-1);
        }

        m_mesh = Mesh::loadFromBinaryStl(filename);

        writeOutput("Loaded triangulation: " + filename);
        writeOutput(QString("Vertices: %1, Edges: %2, Faces: %3").arg(m_mesh->getNumberOfVertices()).arg(m_mesh->getNumberOfEdges()).arg(m_mesh->getNumberOfFaces()));
        if (m_mesh->getNumberOfFaces() == 0) {
            writeOutput("Isosurface Empty!");
            return;
        }

        if (storeImportedMeshAsStl) {
            m_stlWriter->writeMeshFacesToStl(m_mesh, "exact_isosurface_noSmoothing.stl");
            writeOutput(QString("Isosurface stored without smoothing. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    //1 - Do mesh smoothing
    if (m_startPoint <= 1 && m_endPoint >= 1 && doMeshSmoothing) {
        m_mesh->smoothSurface(m_maxSmoothIterations, maxSmoothDistFromOrigin, abortSmoothAverageDist);

        writeOutput(QString("Smoothing done. Time: %1 ms").arg(timeStep.restart()));
        if (storeSmoothedMeshAsStl) {
            m_stlWriter->writeMeshFacesToStl(m_mesh, "exact_isosurface.stl");
            writeOutput(QString("Isosurface stored with smoothing. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    if (m_scaleMesh == 0)    //fit mesh into unit cube
        m_mesh->scaleMeshToFitIntoCube(m_keepMeshProportions, m_scalingFactor);
    else if (m_scaleMesh == 1)  //fit mesh into sphere of diameter 1
        m_mesh->scaleMeshToFitIntoSphere(0.5, QVector3D(0.5, 0.5, 0.5));

    //Export Mesh
    if (exportMesh && ((m_startPoint <= 1 && m_endPoint >= 0) || m_scaleMesh > -1)) {
        m_mesh->saveToFile(meshSaveFilename);
        writeOutput(QString("Mesh exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //2 - Compute the quad mesh
    if (m_startPoint <= 2 && m_endPoint >= 2) {
        QuadmeshGenerator qmGenerator(m_mesh);
        qmGenerator.setParameterizationParameters(param_bt_forQuadmeshGeneration, param_wt, m_parameterInterval, parameterizationStorageReservation_ExpectedEdgesPerVertex);
        if (m_quadmeshAlgorithm == QMA_DelaunayMatching) {
            m_cellMesh = qmGenerator.computeQuadmeshByDelaunayMatching(m_instantTopologyCheck, m_minNumberOfQuadCells, m_minNumberOfInitialCells);
            writeOutput(QString("Quadmesh computed by Delaunay matching! Number of cells: %1. Time: %2 ms").arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else if (m_quadmeshAlgorithm == QMA_VoronoiSubdivision){
            m_cellMesh = qmGenerator.computeQuadmeshBySubdividingVoronoi(m_instantTopologyCheck, m_skipRefinement, m_minNumberOfQuadCells, m_minNumberOfInitialCells);
            writeOutput(QString("Quadmesh computed by Voronoi subdivision! Number of cells: %1. Time: %2 ms").arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else if (m_quadmeshAlgorithm == QMA_CuttingPlanes) {
            QVector<QVector3D> basePoints, normals;

            QString filename = QFileDialog::getOpenFileName(0, QString("Load cutting planes"), m_dmFolder, "*.cp.dm");
            if (filename.isEmpty()) {
                writeOutput(QString("No cutting plane file selected"));
                exit(-1);
            }
            Controller::loadBasepointsAndNormalsFromFile(filename, &basePoints, &normals);

            m_cellMesh = qmGenerator.subdivideMeshByCuttingPlanes(normals, basePoints);
            writeOutput(QString("Quadmesh computed by use of cutting planes, loaded from: %1\n Number of cells: %2. Time: %3 ms").arg(filename).arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else {
            writeOutput("ERROR: No valid quadmesh algorithm chosen!");
            exit(-1);
        }

        //TODO maybe check debugging mode and do regular output if debuging > 1
        qDebug() << "Vertex Valencies" << m_cellMesh->getValencyHistogramCellVertices() << "total:" << m_cellMesh->getNumberOfVertices();
        qDebug() << "Face Valencies" << m_cellMesh->getValencyHistogramCellFaces() << "total:" << m_cellMesh->getNumberOfFaces();

        m_cellMesh->checkAndEnforceCellOrientation();
        m_cellMesh->computeFullParameterization(param_bt_forQuadmeshGeneration);

        if (storeInitialCellMeshAsStl)
            this->cellMeshToStl(QString("initialQuadmesh"), &timeStep, cellMeshIndividualVerticesAndEdgesToStl, cellMeshContentsToStl);
    }

    //3 - Optimize geometry of surface cells
    if (m_startPoint <= 3 && m_endPoint >= 3) {
        for (int i = 0; i < m_numberOfCentralSubdivisions; i++) {
            m_cellMesh->doSubdivisionForFullMesh();
            writeOutput(QString("Cell subdivision %1 / %2 done...").arg(i+1).arg(m_numberOfCentralSubdivisions));
        }

        for (int i = 0; i < m_numberOfVertexOptimizations; i++) {
            m_cellMesh->optimizeAllVertices();
            //m_cellMesh->optimizeAllVerticesSafe();
            writeOutput(QString("Vertex optimization %1 / %2 done...").arg(i+1).arg(m_numberOfVertexOptimizations));
        }

        m_cellMesh->checkAndEnforceCellOrientation();
        writeOutput(QString("Quadmesh optimized. Time: %1 ms").arg(timeStep.restart()));
        if (storeOptimizedCellMeshAsStl)
            this->cellMeshToStl(QString("optimizedSurfaceCellMesh"), &timeStep, cellMeshIndividualVerticesAndEdgesToStl, cellMeshContentsToStl);
    }

//    for (int i = 0; i < m_cellMesh->getNumberOfEdges(); i++)
//        m_cellMesh->getEdge(i)->setG1Edge(false);
//    QVector<int> g1Edges;
//    g1Edges << 21 << 23 << 25 << 47 << 54 << 62 << 58 << 65 << 61 << 33 << 17 << 3 << 49 << 56 << 67 << 68 << 19 << 9 << 1 << 39 << 34 << 32 << 2 << 37 << 16 << 10 << 46 << 40 << 22 << 26 << 31 << 18;
//    g1Edges << 7 << 43 << 41 << 11 << 79 << 75 << 30 << 15 << 76 << 78 << 28 << 77;
//    for (int i = 0; i < g1Edges.size(); i++)
//        m_cellMesh->getEdge(g1Edges[i])->setG1Edge(true);

    //Invert orientation
    if (m_endPoint >= 2 && m_invertCellmeshOrientation) {
        m_cellMesh->invertTotalOrientation();
        writeOutput(QString("Cellmesh orientation inverted"));
    }

    //Export cellmesh
    if (exportCellMesh && (m_startPoint <= 3 || m_invertCellmeshOrientation) && m_endPoint >= 2) {
        m_cellMesh->checkAndEnforceCellOrientation();
        if (QFile::exists(cellMeshSaveFilename)) {
            if (QFile::exists(cellMeshBackupSaveFilename))
                QFile::remove(cellMeshBackupSaveFilename);
            QFile::copy(cellMeshSaveFilename, cellMeshBackupSaveFilename);
        }

        m_cellMesh->saveToFile(cellMeshSaveFilename);
        writeOutput(QString("Cellmesh exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //Export parameterization
    if (exportParameterization && m_startPoint <= 3 && m_endPoint >= 2) {
        if (saveParameterizationsAsImage)
            m_cellMesh->saveParameterizations(parameterizationSaveFilename, parameterizationImageFilename);
        else
            m_cellMesh->saveParameterizations(parameterizationSaveFilename, QString());
        writeOutput(QString("Parameterizations exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //Construct energy matrix
    Eigen::SparseMatrix<double> energyMatrix(m_numberOfControlPoints * m_numberOfControlPoints, m_numberOfControlPoints * m_numberOfControlPoints);
    EnergyMatrixGenerator energyMatGenerator(m_energyFolder, m_loadEnergyMatricesIfPossible, m_storeEnergyMatricesOnHardDrive);
    if (m_startPoint <= 6 && m_endPoint >= 4 && m_weightEnergyTerm > 0.000001) {
        const int numDerivatives = energyDerivatives.size();
        for (int i = 0; i < numDerivatives; i++) {
            const double intervalStart = m_bSplineProperties->getInnerKnotsStart();
            const double intervalEnd = m_bSplineProperties->getInnerKnotsEnd();
            energyMatrix += energyDerivativeFactors[i] * energyMatGenerator.getEnergyMatrixSurface(m_bSplineProperties,
                                                                                                   energyDerivatives[i].first,
                                                                                                   energyDerivatives[i].second,
                                                                                                   intervalStart, intervalEnd, intervalStart, intervalEnd);
        }
        writeOutput(QString("Energy matrix constructed. Time %1 ms").arg(timeStep.restart()));
    }

    //4 - Do B-Spline Approximation
    if (m_startPoint <= 4 && m_endPoint >= 4) {
        //check for consistent orientation
        if (!m_cellMesh->checkForConsistentOrientation())
            writeOutput(QString("WARNING: Individual cell faces are not orientet consistently"));

        //correct parameterizations
        m_cellMesh->computeFullParameterization(param_bt_forSplineFitting);
        const double paramCorrection = m_cellMesh->limitQuadCellParameterizationsToParameterIntervall();
        if (paramCorrection > 0)
            writeOutput(QString("Parameterizations limited to [%1,%2]. Total error: %3").arg(m_parameterInterval.first).arg(m_parameterInterval.second).arg(paramCorrection));

        Splinefitter bSplineFitter(m_cellMesh, m_bSplineProperties);
        bSplineFitter.setEnergyMatrix(energyMatrix);
        bSplineFitter.setWeightEnergyTerm(m_weightEnergyTerm, true);
        bSplineFitter.setNumberOfG1EdgePoints(m_numberOfG1EdgePoints);
        bSplineFitter.setNumberOfG1VertexPoints(m_numberOfG1VertexPoints);
        bSplineFitter.setSizeOfG1VertexDomain(m_g1VertexParameterRange);
        bSplineFitter.setNloMaxIterations(m_nloptMaxIterations);
        bSplineFitter.setNloXPrec(m_nloptXPrec);
        bSplineFitter.setNloFPrec(m_nloptFPrec);
        bSplineFitter.setNloEqualConstPrec(m_nloptConstrPrec);
        bSplineFitter.setNloG0constPrec(m_nloptConstrPrec);
        bSplineFitter.setNloG1constPrec(m_nloptConstrPrec);
        bSplineFitter.setG1VertexConstraintsAsOpt(m_doVertexFittingWithG1Opt);
        bSplineFitter.setG1VertexConstraintsOptFactor(m_vertexFitG1Weight);
        bSplineFitter.setG1EdgeConstraintsAsOpt(m_doEdgeFittingWithG1Opt);
        bSplineFitter.setG1EdgeConstraintsOptFactor(m_edgeFitG1Weight);
        bSplineFitter.setDetailedOutput(m_detailedFittingOutput);

        if (m_useInitialGuess)
            bSplineFitter.setInitialGuess(m_bSplinePatchNetwork);

        //Do B-spline fitting
        if (m_bSplineAlgorithm == BSA_CurveFittingG0) {
            double wEnergy = m_weightEnergyTerm * m_weightEnergyTerm;
            double wLS = 1 - wEnergy;
            m_bSplinePatchNetwork = bSplineFitter.doApproximationLocalC0WithBoundaryInterpolation(true, wLS, wEnergy);
        } else if (m_bSplineAlgorithm == BSA_CurveFittingG0Iterative) {
            QVector<double> wEnergyCurve, wEnergySurf;
            //wEnergyCurve << 0.5 << 0.95 << 0.75 << 0.5 << 0.25 << 0.1 << 0.01 << 0.001 << 0.0001 << 0;
            wEnergyCurve << 0.5 << 0.95 << 0.66 << 0.33 << 0.1 << 0.001 << m_weightEnergyTerm * m_weightEnergyTerm;// << 0.0001 << 0;
            wEnergySurf << 0.5 << 0.95 << 0.75 << 0.5 << 0.25 << 0.1 << 0.01 << m_weightEnergyTerm;// << 0.001 << 0.0001 << 0;
            m_bSplinePatchNetwork = bSplineFitter.doIterativeCurveAndSurfaceFittingC0BoundaryInterpolation(wEnergyCurve, wEnergySurf, true);
        } else if (m_bSplineAlgorithm == BSA_GlobalG0) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationGlobalC0G1Constraints(true);
        } else if (m_bSplineAlgorithm == BSA_SimplifiedGlobalG1) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationGlobalC0G1Constraints(false);
        } else if (m_bSplineAlgorithm == BSA_LocalG1Multiphase) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationLocalG0ConstraintG1Sample();
        } else {
            writeOutput("ERROR: No valid b-spline algorithm chosen!");
            exit(-1);
        }

        writeOutput(QString("B-Spline Approximation done. Algorithm %1. Time: %2 ms").arg(m_bSplineAlgorithm).arg(timeStep.restart()));

        if (storeComputedBSplinesAsStl)
            this->bSplinesToStl(m_stlBSplineSamplePoints, QString("finalBSpline"), &timeStep, bSplineControlPointsToStl, bSplineKnotLinesToStl, bSplineDataErrorToStl, bSplineDataErrorOfIndividualSurfacesToStl);
    }

//    for (int i = 0; i < m_bSplineSurfaces.size(); i++)
//        qDebug() << calculatePatchEnergy(m_bSplineSurfaces[i], energyMatrix, m_numberOfControlPoints);

    //Export B-Splines
    if (exportBSplines && m_startPoint <= 4 && m_endPoint >= 4) {
        if (QFile::exists(bSplineSaveFilename)) {
            if (QFile::exists(bSplineBackupSaveFilename))
                QFile::remove(bSplineBackupSaveFilename);
            QFile::copy(bSplineSaveFilename, bSplineBackupSaveFilename);
        }
        m_bSplinePatchNetwork->saveToFile(bSplineSaveFilename);
        writeOutput(QString("B-Splines exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //5 - Optimize parameters (projection onto tangential space)
    if (m_startPoint <= 5 && m_endPoint >= 5 && doParamOpt) {
        //optimize parameterizations
        for (int iOpt = 0; iOpt < m_numberOfParameterOptimizations; iOpt++) {
            const double improvement = ParameterOptimizer::optimizeParameterizationFullSurfaceNetwork(m_bSplinePatchNetwork, true, m_useNLOforParameterOpt);

            writeOutput(QString("Parameterizations optimized, error new/old: %1 [%2 / %3]. Time: %4 ms").arg(improvement).arg(iOpt + 1).arg(m_numberOfParameterOptimizations).arg(timeStep.restart()));
        }

        //save parameterizations
        if (exportParameterization) {
            if (saveParameterizationsAsImage)
                m_cellMesh->saveParameterizations(parameterizationSaveFilename, parameterizationImageFilename);
            else
                m_cellMesh->saveParameterizations(parameterizationSaveFilename, QString());
            writeOutput(QString("Parameterizations exported. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    //TODO measure:
    //-energy
    //-lenghts of border curves (individual and total)

    //6 - Error measurement
    if (m_startPoint <= 6 && m_endPoint >= 6) {
        //correct parameterizations (in case of new parameterization or step 5 messed something up)
        const double paramCorrection = m_cellMesh->limitQuadCellParameterizationsToParameterIntervall();
        if (paramCorrection > 0)
            writeOutput(QString("Parameterizations limited to [%1,%2]. Total error: %3").arg(m_parameterInterval.first).arg(m_parameterInterval.second).arg(paramCorrection));

        if (m_scaleForEvaluation == 0) {
            double minX, maxX, minY, maxY, minZ, maxZ;
            m_mesh->calculateMinMaxVertexPositions(minX, maxX, minY, maxY, minZ, maxZ);
            m_mesh->scaleMeshToFitIntoCube(m_keepMeshProportions, m_scalingFactor);

            double difX = maxX - minX;
            double difY = maxY - minY;
            double difZ = maxZ - minZ;
            if (m_keepMeshProportions) {
                double max = difX;
                if (difY > max)
                    max = difY;
                if (difZ > max)
                    max = difZ;
                difX = max;
                difY = max;
                difZ = max;
            }

            const QVector3D translationVector(-minX, -minY, -minZ);

            m_cellMesh->translate(translationVector);
            m_cellMesh->scale(m_scalingFactor / difX, m_scalingFactor / difY, m_scalingFactor / difZ);

            m_bSplinePatchNetwork->translate(translationVector);
            m_bSplinePatchNetwork->scale(m_scalingFactor / difX, m_scalingFactor / difY, m_scalingFactor / difZ);
        } else if (m_scaleForEvaluation == 1) {

            double radius;
            QVector3D center;
            m_mesh->calculateMinimumBoundingSphere(radius, center);

            const double scalingFactor = m_scalingFactor/(2 * radius);
            const QVector3D targetCenter(m_scalingFactor/2, m_scalingFactor/2, m_scalingFactor/2);
            m_mesh->scaleMeshToFitIntoSphere(m_scalingFactor/2, targetCenter);

            m_cellMesh->translate(-center);
            m_cellMesh->scale(scalingFactor, scalingFactor, scalingFactor);
            m_cellMesh->translate(targetCenter);

            m_bSplinePatchNetwork->translate(-center);
            m_bSplinePatchNetwork->scale(scalingFactor, scalingFactor, scalingFactor);
            m_bSplinePatchNetwork->translate(targetCenter);

        }

        //evaluate
        SplineEvaluator evaluator(m_bSplinePatchNetwork, &energyMatrix, numberOfEvalPointsForG1, doOptimizationForDataError);

        evaluator.writeQualityFile(qualityFilename);

        writeOutput(QString("Spline area: %1 Energy: %2 Energy/Area: %3").arg(evaluator.getTotalSurfaceAreaBSplines()).arg(evaluator.getTotalEnergy()).arg(evaluator.getTotalEnergy()/evaluator.getTotalSurfaceAreaBSplines()));
        writeOutput(QString("RMS Error: %1 Max Error: %2").arg(evaluator.getRmsError()).arg(evaluator.getMaxDataError()));
        writeOutput(QString("BBNormalized: %1 BBNMax: %2").arg(evaluator.getBoundingBoxScaledRmsError()).arg(evaluator.getBoundingBoxScaledMaxDataError()));
        writeOutput(QString("BSNormalized [%]: %1 BSNMax: %2").arg(evaluator.getBoundingSphereScaledRmsError()).arg(evaluator.getBoundingSphereScaledMaxDataError()));
        writeOutput(QString("Average G1 Angle Error: %1 Maximum Error: %2").arg(evaluator.getTotalAverageAngleError()).arg(evaluator.getMaxAngleError()));
        writeOutput(QString("B-Spline quality measured. Time: %1 ms").arg(timeStep.restart()));

        //matlab export
        if (m_doMatlabOutput) {
            Controller::removeRecursively(m_dmFolder + matlabSubfolder);
            if (!QDir(m_dmFolder + matlabSubfolder).exists())
                QDir().mkpath(m_dmFolder + matlabSubfolder);

            evaluator.writeMatlabTriangulationStructureFile(matlabTriStructureFilename);
            evaluator.writeMatlabDataErrorFiles(matlabTriErrorFilename, matlabBSErrorFilename);
            evaluator.writeMatlabG1ErrorAtPLVFiles(matlabTriBordersFilename, matlabBSBorderFilename);
            evaluator.writeMatlabRegularSamplingStructureFile(numberOfEvalPointsForCurvature, matlabRegularSamplingStructureFilename);
            evaluator.writeMatlabRegularSampledCurvature(numberOfEvalPointsForCurvature, matlabRegularSamplingCurvatureFilename);
            evaluator.writeMatlabRegularSamplingG1Error(matlabRegularSamplingBorderFilename);

            writeOutput(QString("Matlab export done: %1 ms").arg(timeStep.restart()));
        }

        if (m_do3DImageOutput) {
            int imageSizeX = m_3DimageSizeX;
            int imageSizeY = m_3DimageSizeY;
            int imageSizeZ = m_3DimageSizeZ;
            QVector3D translationVector(0, 0, 0);
            if (m_3DimageSizeX == -1) {
                int minX, maxX, minY, maxY, minZ, maxZ;
                m_bSplinePatchNetwork->getBoundingBox(minX, maxX, minY, maxY, minZ, maxZ);
                //translate by -(minX, minY, minZ), so that lower corner of the bounding box is (0, 0, 0)
                //add padding in each direction
                translationVector = QVector3D(iassPadding - minX, iassPadding - minY, iassPadding - minZ);
                imageSizeX = 2 * iassPadding + maxX - minX;
                imageSizeY = 2 * iassPadding + maxY - minY;
                imageSizeZ = 2 * iassPadding + maxZ - minZ;
            }

            m_bSplinePatchNetwork->translate(translationVector);
            Voxelizer::splinesToVoxelImage(m_bSplinePatchNetwork, imageSizeX, imageSizeY, imageSizeZ, iassOutputFilename, iassSurfaceSamplePoints);
            m_bSplinePatchNetwork->translate(-translationVector);

            writeOutput(QString("Volume image done. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    if (m_startPoint <= 7 && m_endPoint >= 7) {
        //Generation of view point candidates (VPC) for surface inspection

        //const int m_mbsiMetricType = 1; //TODO replace with enum later and make it a parameter
        const int maxSubdivisionDepth = 50; //TODO make that a parameter

        const int precomputeDepth = 1;  //Only used by thin plate metric
        const int surfaceMetricNumericalSampleRatePerInterval = 2 * m_splineOrder;
        const bool integrateCurvatureIn3D = true;

        QString factorToken("a");
        if (m_mbsiThresholdIsAVGFactor)
            factorToken = "f";
        if (m_mbsiTargetNumberOfVPC > -1) {
            factorToken = "n";
            m_mbsiThreshold = m_mbsiTargetNumberOfVPC;
        }

        //Set up metric and compute average for entire surface network
        SubdivisionAbstractMetric *metric = 0;
        if (m_mbsiMetricType == 0) {
            metric = new SubdivisionThinPlateEnergyMetric(&energyMatGenerator, m_bSplineProperties);
            ((SubdivisionThinPlateEnergyMetric *) metric)->precomputeThinPlateMatrices(precomputeDepth);
        } else if (m_mbsiMetricType == 1) {
            metric = new SubdivisionExactCurvatureMetric(surfaceMetricNumericalSampleRatePerInterval, integrateCurvatureIn3D);
        } else if (m_mbsiMetricType == 2) {
            metric = new SubdivisionSurfaceAreaMetric(surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 3) {
            metric = new SubdivisionNormalDeviationMetric(surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 4) {
            metric = new SubdivisionInspectionAngleMetric(m_mbsiCamDistance, surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 5) {
            //Stretched cylinder
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(21, -9, 18), 5, 1);
        } else if (m_mbsiMetricType == 6) {
            //Spring
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(4, 4, 20), 3, 1);
        } else if (m_mbsiMetricType == 7) {
            //Tangle cube
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(70, 70, 70), 20, 1);
        } else {
            qDebug() << "MBSI: No valid metric specified!";
            exit(-1);
        }
//            srand(time(NULL));
//            const int numberOfRandomHighlights = 5;
//            const double randomHighlighRadius = 0.15;
//            for (int i = 0; i < numberOfRandomHighlights; i++) {
//                const int patchId = rand() % m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
//                const double randomU = m_parameterInterval.first + (((double) rand() / (RAND_MAX)) * (m_parameterInterval.second - m_parameterInterval.first));
//                const double randomV = m_parameterInterval.first + (((double) rand() / (RAND_MAX)) * (m_parameterInterval.second - m_parameterInterval.first));
//                const QVector3D point = m_bSplinePatchNetwork->getBSplineSurface(patchId)->evaluate(randomU, randomV);

//                qDebug() << patchId << randomU << randomV;
//                qDebug() << point << randomHighlighRadius;

//                ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(point, randomHighlighRadius, 1);
//            }
        const QString metricNameToken = metric->getMetricNameToken();

        //Define threshold
        const int numberOfSurfaces = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
        double averageMetric = 1;
        if (m_mbsiThresholdIsAVGFactor && m_mbsiTargetNumberOfVPC < 0) {
            averageMetric = 0;
            //TODO these values get recomputed at depth = 0 for the subdivision
            for (int i = 0; i < numberOfSurfaces; i++) {
                const double value = metric->evaluate(m_bSplinePatchNetwork->getBSplineSurface(i),
                                                      m_parameterInterval.first, m_parameterInterval.second,
                                                      m_parameterInterval.first, m_parameterInterval.second);
                averageMetric += value;
                qDebug() << i << value;
            }
            averageMetric /= (double) numberOfSurfaces;
            qDebug() << "Total:" << averageMetric * numberOfSurfaces << "Average:" << averageMetric;
        }
        const double subdivisionThreshold = m_mbsiThreshold * averageMetric;

        writeOutput(QString("Metric and threshold setup done. Absolute threshold value: %1, time: %2").arg(subdivisionThreshold).arg(timeStep.restart()));


        //Do the subdivision
        QTime timeSubdivision;
        timeSubdivision.start();

        QVector<SubSurfaceTree *> subdivisions;
        if (m_mbsiTargetNumberOfVPC > -1)
            subdivisions = SurfaceSubdivider::subdivideSurfacesUntilFixedNumber(m_bSplinePatchNetwork->getBSplineSurfaces(), metric, 100, maxSubdivisionDepth);
        else
            subdivisions = SurfaceSubdivider::subdivideSurfacesThresholdBased(m_bSplinePatchNetwork->getBSplineSurfaces(), metric, subdivisionThreshold, maxSubdivisionDepth);

        writeOutput(QString("Subdivision done! Time %1 ms").arg(timeSubdivision.restart()));

        //Write subdivisions to stl
        if (writeMBSISubdivisionsToSTL)
            m_stlWriter->writeBSplineSubdivisionsToStl(m_bSplinePatchNetwork->getBSplineSurfaces(), subdivisions, QString("subdivisions_%1_th%2%3_LW%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_stlLineWidth), m_stlLineWidth, m_stlBSplineSamplePoints);

        //Set up vectors for output
        QVector<QVector3D> pivotPointsCurv;
        QVector<QVector3D> camPositionsCurv;
        QVector<int> depths;

        //Calculate VPC
        ViewpointCandidateGenerator viewpointGenerator(m_bSplinePatchNetwork);
        viewpointGenerator.generateViewPoints_subdivision(camPositionsCurv, pivotPointsCurv, depths, m_mbsiCamDistance, subdivisions);

        //Process result
        const int numSamplePoints = camPositionsCurv.size();
        QVector<QVector3D> normalsCurv(numSamplePoints);
        QVector<Edge *> camDirectionsCurv(numSamplePoints);
        for (int i = 0; i < numSamplePoints; i++) {
            normalsCurv[i] = (camPositionsCurv[i] - pivotPointsCurv[i]).normalized();   //This is not computed by the VPC generator anymore
            Vertex *v0 = new Vertex(pivotPointsCurv[i], -1);
            Vertex *v1 = new Vertex(pivotPointsCurv[i] + m_mbsiCamDirVectorLength * normalsCurv[i], -1);
            camDirectionsCurv[i] = new Edge(v0, v1, -1);
        }

        if (writeMBSIViewpointCandidatesToSTL) {
            m_stlWriter->writePoints(&camPositionsCurv, m_mbsiCamPointSize, QString("vpc_cam_%1_th%2%3_size%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsCurv, m_mbsiPivotPointSize, QString("vpc_piv_%1_th%2%3_size%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsCurv, QString("vpc_normals_%1_th%2%3_size%4_dist%5.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiCamDirVectorSize).arg(m_mbsiCamDirVectorLength), m_mbsiCamDirVectorSize);
        }

        std::vector<std::vector<double> > positionList(numSamplePoints);
        for (int i = 0; i < numSamplePoints; i++) {
            std::vector<double> point(10);
            //Cam position
            point[0] = camPositionsCurv[i].x();
            point[1] = camPositionsCurv[i].y();
            point[2] = camPositionsCurv[i].z();
            //Cam direction is the negative of the normal vector
            const QVector3D normal = normalsCurv[i];
            point[3] = -normal.x();
            point[4] = -normal.y();
            point[5] = -normal.z();
            //Orientation (up vector)
            QVector3D up = QVector3D(0, 0, 1) - QVector3D::dotProduct(QVector3D(0, 0, 1), normal) * normal;
            if (up == QVector3D(0, 0, 0))
                up = QVector3D(0, 1, 0);
            up.normalize();
            point[6] = up.x();
            point[7] = up.y();
            point[8] = up.z();

            //depth of the subdivision
            point[9] = depths[i];

            positionList[i] = point;
        }

        QString positionFilename(m_dmFolder + QString("VPC_%1_dist%2_th%3%4.txt").arg(metricNameToken).arg(m_mbsiCamDistance).arg(m_mbsiThreshold).arg(factorToken));
        std::string stdFilename = positionFilename.toStdString();
        writePositions(positionList, stdFilename);

        writeOutput(QString("Subdivision based sample points done. Metric: %1, number of points %2. Time: %3").arg(metricNameToken).arg(numSamplePoints).arg(timeStep.restart()));

        if ((m_mbsiMetricType == 0 || m_mbsiMetricType == 1 || m_mbsiMetricType == 4) && doExtendedMBSIOutput) {
            QTime t;
            t.start();
            qDebug() << "starting timing for new plv output";
            BSplineTriangulation *triangulation = BSplineTriangulation::createRegularSampledTriangulation(m_bSplinePatchNetwork->getNumberOfBSplineSurfaces(), m_stlBSplineSamplePoints,
                                                                                                          m_parameterInterval.first, m_parameterInterval.second, m_parameterInterval.first, m_parameterInterval.second,
                                                                                                          true);    //TODO test with false as well
            qDebug() << "triangulation done" << t.restart();
            ColorMapGenerator colMapGenerator(true);    //TODO test with different fixed color ranges
            colMapGenerator.setColorMapToHeat(false);
            QVector<QColor> colorMap = colMapGenerator.createColorMap(m_bSplinePatchNetwork->getBSplineSurfaces(), triangulation, metric);
            qDebug() << "color map done" << t.restart();
            PlyWriter::writeBSplineTriangulationWithVertexColors(m_bSplinePatchNetwork, triangulation, colorMap, m_dmFolder + QString("metricValues_%1.ply").arg(metricNameToken));
            qDebug() << "ply output done" << t.restart();
            QString matlabPlotCommand = colMapGenerator.generateMatlabPlotCommand(101, QString("./cbar_%1.png").arg(metricNameToken));
            std::cout << matlabPlotCommand.toStdString() << std::endl;
            QFile file(m_dmFolder + QString("colorbar_%1.m").arg(metricNameToken));
            file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
            QTextStream out(&file);
            out << matlabPlotCommand;
            file.close();
            qDebug() << "matlab colorbar done" << t.restart();
    //        for (int i = 0; i < numberOfSurfaces; i++) {
    //            qDebug() << "surface" << i;
    //            this->createMatlab2DOutputForMetric(m_bSplinePatchNetwork->getBSplineSurface(i), QString("../out/metric_%1_%2.csv").arg(metricNameToken).arg(i), metric, subdivisions[i], surfaceMetricNumericalSampleRate, m_parameterInterval);
    //        }
    //        qDebug() << "2D matlab output done" << t.restart();
        }

        delete metric;
        foreach (SubSurfaceTree *subsurface, subdivisions)
            delete subsurface;

        //Sample point generation with reference algorithms
        if (m_mbsiDoReferenceAlgorithms) {
            QVector<QVector3D> camPositionsVertexBased, pivotPointsVertexBased, camPositionsRandomCenter, pivotPointsRandomCenter, camPositionsRandomVertex, pivotPointsRandomVertex;
            viewpointGenerator.generateViewPoints_useVertices(camPositionsVertexBased, pivotPointsVertexBased, m_mbsiCamDistance);
            viewpointGenerator.generateViewPoints_random(camPositionsRandomCenter, pivotPointsRandomCenter, m_mbsiCamDistance, numSamplePoints, false);
            viewpointGenerator.generateViewPoints_random(camPositionsRandomVertex, pivotPointsRandomVertex, m_mbsiCamDistance, numSamplePoints, true);

            QVector<Edge *> camDirectionsVertexBased(camPositionsVertexBased.size());
            for (int i = 0; i < camPositionsVertexBased.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsVertexBased[i], -1);
                Vertex *v1 = new Vertex(pivotPointsVertexBased[i], -1);
                camDirectionsVertexBased[i] = new Edge(v0, v1, -1);
            }
            QVector<Edge *> camDirectionsRandomCenter(camPositionsRandomCenter.size());
            for (int i = 0; i < camPositionsRandomCenter.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsRandomCenter[i], -1);
                Vertex *v1 = new Vertex(pivotPointsRandomCenter[i], -1);
                camDirectionsRandomCenter[i] = new Edge(v0, v1, -1);
            }
            QVector<Edge *> camDirectionsRandomVertex(camPositionsRandomVertex.size());
            for (int i = 0; i < camPositionsRandomVertex.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsRandomVertex[i], -1);
                Vertex *v1 = new Vertex(pivotPointsRandomVertex[i], -1);
                camDirectionsRandomVertex[i] = new Edge(v0, v1, -1);
            }
            m_stlWriter->writePoints(&camPositionsVertexBased, m_mbsiCamPointSize, QString("mbsi_vertex_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsVertexBased, m_mbsiPivotPointSize, QString("mbsi_vertex_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsVertexBased, QString("mbsi_vertex_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            m_stlWriter->writePoints(&camPositionsRandomCenter, m_mbsiCamPointSize, QString("mbsi_rndCent_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsRandomCenter, m_mbsiPivotPointSize, QString("mbsi_rndCent_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsRandomCenter, QString("mbsi_rndCent_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            m_stlWriter->writePoints(&camPositionsRandomVertex, m_mbsiCamPointSize, QString("mbsi_rndVert_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsRandomVertex, m_mbsiPivotPointSize, QString("mbsi_rndVert_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsRandomVertex, QString("mbsi_rndVert_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            writeOutput(QString("Sample points based on reference algorithms done. Time: %2").arg(timeStep.restart()));
        }

    }
    writeOutput(QString("All done! Total Time: %1 ms").arg(timeTotal.elapsed()));
}


void Controller::setStlWriter(StlWriter *writer)
{
    m_stlWriter = writer;
}

void Controller::setDebugingMode(int mode)
{
    m_debugingMode = mode;
}

void Controller::setConsoleOutput(bool mode)
{
    m_consoleOutput = mode;
}

void Controller::setLogOutput(bool mode)
{
    m_logOutput = mode;
}

void Controller::setStartPoint(int startPoint)
{
    m_startPoint = startPoint;
}

void Controller::setEndPoint(int endPoint)
{
    m_endPoint = endPoint;
}

void Controller::setDmFolder(QString dmFolder)
{
    m_dmFolder = dmFolder;
    m_logFilename = QString(m_dmFolder + "log.txt.dm");
}

void Controller::setEnergyFolder(QString energyFolder)
{
    m_energyFolder = energyFolder;
}

void Controller::setStlBSplineSamplePoints(int samplePoints)
{
    m_stlBSplineSamplePoints = samplePoints;
}

void Controller::setStlLineWidth(double lineWidth)
{
    m_stlLineWidth = lineWidth;
}

void Controller::setNumberOfThreads(int numberOfThreads)
{
    m_fixedNumberOfThreads = numberOfThreads;
}

void Controller::setParameterInterval(double a, double b)
{
    m_parameterInterval = QPair<double, double>(a, b);
    m_useLoadedParamInterval = false;
}

void Controller::setMeshScaling(int scaling)
{
    m_scaleMesh = scaling;
}

void Controller::setKeepingProportions(bool keepProportions)
{
    m_keepMeshProportions = keepProportions;
}

void Controller::setScalingFactor(double scalingFactor)
{
    m_scalingFactor = scalingFactor;
}

void Controller::setmaxSmoothIterations(int maxIterations)
{
    m_maxSmoothIterations = maxIterations;
}

void Controller::setQuadmeshAlgorithm(Controller::QuadmeshAlgorithm alg)
{
    m_quadmeshAlgorithm = alg;
}

void Controller::setMinNumberOfInitialCells(int numberOfCells)
{
    m_minNumberOfInitialCells = numberOfCells;
}

void Controller::setMinNumberOfQuadCells(int numberOfCells)
{
    m_minNumberOfQuadCells = numberOfCells;
}

void Controller::setInstantTopologyCheck(bool instantTopologyCheck)
{
    m_instantTopologyCheck = instantTopologyCheck;
}

void Controller::setSkipRefinement(bool skipRefinement)
{
    m_skipRefinement = skipRefinement;
}

void Controller::setNumberOfCentralSubdivisions(int numberOfSubdivisions)
{
    m_numberOfCentralSubdivisions = numberOfSubdivisions;
}

void Controller::setNumberOfVertexOptimizations(int numberOfOptimizations)
{
    m_numberOfVertexOptimizations = numberOfOptimizations;
}

void Controller::setBSplineAlgorithm(Controller::BSplineAlgorithm alg)
{
    m_bSplineAlgorithm = alg;
}

void Controller::setRecomputeParameterization(bool recomputeParam)
{
    m_recomputeParameterization = recomputeParam;
}

void Controller::setInvertCellmeshOrientation(bool invert)
{
    m_invertCellmeshOrientation = invert;
}

void Controller::setUsageOfInitialGuess(bool useInitialGuess)
{
    m_useInitialGuess = useInitialGuess;
}

void Controller::setNumberOfControlPoints(int numberOfControlPointsPerRow)
{
    m_numberOfControlPoints = numberOfControlPointsPerRow;
}

void Controller::setSplineOrder(int splineOrder)
{
    m_splineOrder = splineOrder;
}

void Controller::setWeightEnergyTerm(double weight)
{
    m_weightEnergyTerm = weight;
}

void Controller::setG1VertexParameterRange(double range)
{
    m_g1VertexParameterRange = range;
    m_g1VertexParameterRangeSet = true;
}

void Controller::setNumberOfG1VertexPoints(int numberOfPoints)
{
    m_numberOfG1VertexPoints = numberOfPoints;
}

void Controller::setVertexFittingG1Opt(bool mode)
{
    m_doVertexFittingWithG1Opt = mode;
}

void Controller::setVertexFitG1Weight(double weight)
{
    m_vertexFitG1Weight = weight;
}

void Controller::setNumberOfG1EdgePoints(int numberOfPoints)
{
    m_numberOfG1EdgePoints = numberOfPoints;
}

void Controller::setEdgeFittingG1Opt(bool mode)
{
    m_doEdgeFittingWithG1Opt = mode;
}

void Controller::setEdgeFitG1Weight(double weight)
{
    m_edgeFitG1Weight = weight;
}

void Controller::setDetailedFittingOutput(bool detailedFittingOutput)
{
    m_detailedFittingOutput = detailedFittingOutput;
}

void Controller::setEnergyMatrixStorageOnHarddrive(bool storeMatrices)
{
    m_storeEnergyMatricesOnHardDrive = storeMatrices;
}

void Controller::setEnergyMatrixLoadingFromHarddrive(bool loadMatrices)
{
    m_loadEnergyMatricesIfPossible = loadMatrices;
}

void Controller::setNumberOfParameterOptimizations(int numberOfOptimizations)
{
    m_numberOfParameterOptimizations = numberOfOptimizations;
}

void Controller::setUseOfNLOforParamOpt(bool useNLO)
{
    m_useNLOforParameterOpt = useNLO;
}

void Controller::setScalingForEvaluation(int scaling)
{
    m_scaleForEvaluation = scaling;
}

void Controller::setMatlabOutput(bool doMatlabOutput)
{
    m_doMatlabOutput = doMatlabOutput;
}

void Controller::set3DImageOutput(bool do3DImage)
{
    m_do3DImageOutput = do3DImage;
}

void Controller::set3DImageSizes(int sizeX, int sizeY, int sizeZ)
{
    m_3DimageSizeX = sizeX;
    m_3DimageSizeY = sizeY;
    m_3DimageSizeZ = sizeZ;
}

void Controller::setMbsiMetricType(int metricType)
{
    m_mbsiMetricType = metricType;
}

void Controller::setMbsiThreshold(double thresholdFactor, bool factorOfAverage)
{
    m_mbsiThreshold = thresholdFactor;
    m_mbsiThresholdIsAVGFactor = factorOfAverage;
}

void Controller::setMbsiTargetNumberOfVPC(int number)
{
    m_mbsiTargetNumberOfVPC = number;
}

void Controller::setMbsiDoReferenceAlgorithms(bool doAlgorithms)
{
    m_mbsiDoReferenceAlgorithms = doAlgorithms;
}

void Controller::setMbsiPivotPointSize(double size)
{
    m_mbsiPivotPointSize = size;
}

void Controller::setMbsiCamPointSize(double size)
{
    m_mbsiCamPointSize = size;
}

void Controller::setMbsiCamDirectionSize(double size)
{
    m_mbsiCamDirVectorSize = size;
}

void Controller::setMbsiCamDirectionLength(double length)
{
    m_mbsiCamDirVectorLength = length;
}

void Controller::setMbsiCamDistance(double distance)
{
    m_mbsiCamDistance = distance;
}

void Controller::setNloptMaxIterations(int maxIterations)
{
    m_nloptMaxIterations = maxIterations;
}

void Controller::setNloptFPrec(double fPrec)
{
    m_nloptFPrec = fPrec;
}

void Controller::setNloptXPrec(double xPrec)
{
    m_nloptXPrec = xPrec;
}

void Controller::setNloptConstrPrec(double constrPrec)
{
    m_nloptConstrPrec = constrPrec;
}

bool Controller::loadBasepointsAndNormalsFromFile(QString filename, QVector<QVector3D> *basepoints, QVector<QVector3D> *normals)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly))
        return false;
    QTextStream in(&file);

    QVector<QVector3D> *current = basepoints;
    while (!in.atEnd()) {
        QString line = in.readLine();

        if (line == "===") {
            current = normals;
        } else {
            QStringList elements = line.split(QString(" "), QString::SkipEmptyParts);
            if (elements.size() != 3)
                return false;

            current->append(QVector3D(elements[0].toFloat(), elements[1].toFloat(), elements[2].toFloat()));
        }
    }
    if (basepoints->size() == normals->size())
        return true;
    else
        return false;
}

void Controller::createMatlab2DOutputForMetric(BSplineSurface *surf, QString csvFilename, SubdivisionAbstractMetric *metric, SubSurfaceTree *subdivision, const int numberOfEvalPoints, QPair<double, double> parameterInterval)
{
    QFile file(csvFilename);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);

    const double uDiff = parameterInterval.second - parameterInterval.first;
    const double vDiff = parameterInterval.second - parameterInterval.first;
    for (int j = 0; j < numberOfEvalPoints; j++) {
        const double v = parameterInterval.first + vDiff * (double) j / (numberOfEvalPoints - 1);
        for (int i = 0; i < numberOfEvalPoints; i++) {
            const double u = parameterInterval.first + uDiff * (double) i / (numberOfEvalPoints - 1);
            const double value = metric->evaluatePoint(surf, u, v);
            //out << u << ", " << v << ", " << value << "\n:";
            out << value;
            if (i == numberOfEvalPoints - 1)
                out << "\n";
            else
                out << ", ";
        }
    }

    qDebug() << QString("[U V] = meshgrid(1:%1:1,0:%1:1);").arg((double) 1/ (numberOfEvalPoints - 1));
    //pcolor(U, V, C)
    //shading flat

    file.close();

    QVector<SubSurfaceTree *> leaves = subdivision->getAllLeaves();
    const int numLeaves = leaves.size();
    const double uMax = subdivision->getU1();
    const double vMax = subdivision->getV1();

    for (int iLeaf = 0; iLeaf < numLeaves; iLeaf++) {
        SubSurfaceTree *leaf = leaves[iLeaf];
        const double u0 = leaf->getU0();
        const double u1 = leaf->getU1();
        const double v0 = leaf->getV0();
        const double v1 = leaf->getV1();

        //line (u0, v1) to (u1, v1)
        if (v1 != vMax)
            qDebug() << "plot([" << u0 << u1 << "],[" << v1 << v1 << "], 'b')";

        //line (u1, v0) to (u1, v1)
        if (u1 != uMax)
            qDebug() << "plot([" << u1 << u1 << "],[" << v0 << v1 << "], 'b')";
    }
}

//void SplineGenerator::createDebugSpline()
//{
//    QVector<double> knots;
//    knots << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1;
//    QVector<QVector<QVector3D> > controlPoints(4, QVector<QVector3D>(4));
//    controlPoints[0][0] = QVector3D(00, 00, 00);
//    controlPoints[0][1] = QVector3D(00, 10, 00);
//    controlPoints[0][2] = QVector3D(00, 20, 00);
//    controlPoints[0][3] = QVector3D(00, 30, 00);

//    controlPoints[1][0] = QVector3D(10, 00, 00);
//    controlPoints[1][1] = QVector3D(10, 10, 10);
//    controlPoints[1][2] = QVector3D(10, 20, 10);
//    controlPoints[1][3] = QVector3D(10, 30, 00);

//    controlPoints[2][0] = QVector3D(20, 00, 00);
//    controlPoints[2][1] = QVector3D(20, 10, 20);
//    controlPoints[2][2] = QVector3D(20, 20, 20);
//    controlPoints[2][3] = QVector3D(20, 30, 00);

//    controlPoints[3][0] = QVector3D(30, 00, 00);
//    controlPoints[3][1] = QVector3D(30, 10, 20);
//    controlPoints[3][2] = QVector3D(30, 20, 20);
//    controlPoints[3][3] = QVector3D(30, 30, 00);

//    BSplineSurface surf;
//    surf.setKnotsU(knots);
//    surf.setKnotsV(knots);
//    surf.setControlPoints(controlPoints);
//    surf.setOrder(3);
//    m_bSplineSurfaces.append(&surf);

//    controlPoints = QVector<QVector<QVector3D> >(4, QVector<QVector3D>(4));
//    controlPoints[0][0] = QVector3D(30, 00, 00);
//    controlPoints[0][1] = QVector3D(30, 10, 20);
//    controlPoints[0][2] = QVector3D(30, 20, 20);
//    controlPoints[0][3] = QVector3D(30, 30, 00);

//    controlPoints[1][0] = QVector3D(40, 00, 00);
//    controlPoints[1][1] = QVector3D(40, 10, 20);
//    controlPoints[1][2] = QVector3D(40, 20, 20);
//    controlPoints[1][3] = QVector3D(40, 30, 00);

//    controlPoints[2][0] = QVector3D(50, 00, 00);
//    controlPoints[2][1] = QVector3D(50, 10, 10);
//    controlPoints[2][2] = QVector3D(50, 20, 10);
//    controlPoints[2][3] = QVector3D(50, 30, 00);

//    controlPoints[3][0] = QVector3D(60, 00, 00);
//    controlPoints[3][1] = QVector3D(60, 10, 00);
//    controlPoints[3][2] = QVector3D(60, 20, 00);
//    controlPoints[3][3] = QVector3D(60, 30, 00);
//    BSplineSurface surf2;
//    surf2.setKnotsU(knots);
//    surf2.setKnotsV(knots);
//    surf2.setControlPoints(controlPoints);
//    surf2.setOrder(3);
//    m_bSplineSurfaces.append(&surf2);

//    BSplineCurve curv;
//    curv.setKnots(knots);
//    curv.setOrder(3);
//    curv.setControlPoints(controlPoints[0]);
//    m_bSplineCurves.append(&curv);

//    this->bSplinesToStl(51, QString("demo"), &time);
//}

void Controller::writeOutput(const QString msg) const
{
    writeToLogOnly(msg);
    writeToConsoleOnly(msg);
}

void Controller::writeToConsoleOnly(const QString msg) const
{
    if (m_consoleOutput)
        std::cout << msg.toStdString() << std::endl;
}

void Controller::writeToLogOnly(const QString msg) const
{
    if (m_logOutput) {
        if (!QDir(m_dmFolder).exists())
            QDir().mkpath(m_dmFolder);
        QFile logFile(m_logFilename);
        if (logFile.exists())
            logFile.open(QIODevice::Append | QIODevice::Text);
        else
            logFile.open(QIODevice::WriteOnly | QIODevice::Text);

        QTextStream out(&logFile);

        out << msg << "\n";

        logFile.close();
    }
}

void Controller::clearLogFile()
{
    QFile logFile(m_logFilename);
    logFile.open(QFile::WriteOnly | QFile::Truncate);
    logFile.close();
}

bool Controller::removeRecursively(const QString &dirName)
{
    //QDir offers this function since QT 5.0
    bool result = true;
    QDir dir(dirName);

    if (dir.exists()) {
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
            if (info.isDir()) {
                result = removeRecursively(info.absoluteFilePath());
            } else {
                result = QFile::remove(info.absoluteFilePath());
            }

            if (!result) {
                return result;
            }
        }
        result = dir.rmdir(dirName);
    }
    return result;
}

void Controller::cellMeshToStl(QString basename, QTime *timeStep, bool storeIndividualCellVerticesAndEdges, bool storeCellContents)
{
    m_cellMesh->cleanUpMeshIndices();   //Neccessary to avoid NULLs in vertex, edge, and face lists

    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    QVector<Edge *> surfaceCellMeshEdges;
    for (int i = 0; i < numberOfCellEdges; i++) {
        if (m_cellMesh->getEdge(i)) {
            QVector<PolylineVertex *> *polyVertices = m_cellMesh->getEdge(i)->getPolylineVertices();
            for (int j = 1; j < polyVertices->size(); j++) {
                surfaceCellMeshEdges.append(new Edge(polyVertices->at(j-1), polyVertices->at(j), -1));
            }
        }
    }
    m_stlWriter->writeEdgesToStl(&surfaceCellMeshEdges, QString(basename).append(QString("_borders_LW%1.stl").arg(m_stlLineWidth)), m_stlLineWidth);
    foreach (Edge *e, surfaceCellMeshEdges)
        delete e;
    writeOutput(QString(basename).append(QString(" - borders stored. Time: %1 ms").arg(timeStep->restart())));

    if (storeIndividualCellVerticesAndEdges) {
        for (int i = 0; i < m_cellMesh->getNumberOfVertices(); i++) {
            QVector<QVector3D> tmp;
            tmp << m_cellMesh->getVertex(i)->getPosition();
            m_stlWriter->writePoints(&tmp, 0.5, QString("debug/").append(basename).append(QString("_CellVertex_%1.stl").arg(i)));
        }
        writeOutput(QString(basename).append(QString(" - cell vertices stored. Time: %1 ms").arg(timeStep->restart())));
        for (int i = 0; i < m_cellMesh->getNumberOfEdges(); i++) {
            CellEdge *ce = m_cellMesh->getEdge(i);
            QVector<Edge *> tmp;
            for (int j = 1; j < ce->getPolylineVertices()->size(); j++)
                tmp << new Edge(ce->getPolylineVertices()->at(j - 1), ce->getPolylineVertices()->at(j), -1);
            m_stlWriter->writeEdgesToStl(&tmp, QString("debug/").append(basename).append(QString("_CellEdge_%1.stl").arg(i)), m_stlLineWidth);
            for (int i = 0; i < tmp.size(); i++)
                delete tmp[i];
        }
        writeOutput(QString(basename).append(QString(" - cell edges stored. Time: %1 ms").arg(timeStep->restart())));
    }

    if (storeCellContents) {
        for (int i = 0; i < m_cellMesh->getNumberOfFaces(); i++) {
            if (m_cellMesh->getFace(i)) {
                m_stlWriter->writeCellFaceContentToStl(m_cellMesh->getFace(i), QString("debug/").append(basename).append(QString("_CellFace_Content_%1.stl").arg(i)), 0.2, 0.1);
                m_stlWriter->writeCellBorderContentToStl(m_cellMesh->getFace(i), m_mesh, QString("debug/").append(basename).append(QString("_CellFace_Border_%1.stl").arg(i)), 0.2, 0.1);
                m_cellMesh->getCellFaceParameterization(i)->saveAsImage(QString(m_stlWriter->getDefaultFolder()).append(QString("debug/")).append(basename).append(QString("_Param_%1.bmp").arg(i)), 1000, 2, 4);
            }
        }
        writeOutput(QString(basename).append(QString(" - individual cell content stored. Time: %1 ms").arg(timeStep->restart())));
    }
}

void Controller::bSplinesToStl(int numberOfEvalPoints, QString basename, QTime *timeStep, bool doControlPoints, bool doKnotLines, bool doDataError, bool doDataErrorForIndividualSurfaces)
{
    m_stlWriter->writeBSplineCurvesToStl(m_bSplinePatchNetwork->getBSplineCurves(), QString(basename).append("_CurveNetwork_LW%1.stl").arg(m_stlLineWidth), m_stlLineWidth, numberOfEvalPoints);
    m_stlWriter->writeBSplineSurfacesToStl(m_bSplinePatchNetwork->getBSplineSurfaces(), QString(basename).append(QString("_SurfaceNetwork.stl")), numberOfEvalPoints);

    writeOutput(QString(basename).append(QString(" - surface patch network stored. Time: %1 ms").arg(timeStep->restart())));

    if (doControlPoints) {
        m_stlWriter->writeBSplineControlPointsToStl(m_bSplinePatchNetwork, QString(basename).append("_controlPointGrid.stl"), QString(basename).append("controlPoints_LW%1.stl").arg(m_stlLineWidth), m_stlLineWidth, m_stlLineWidth);
        m_stlWriter->writeControlPointsOfIndividualPhases(m_bSplinePatchNetwork, basename + "_phaseCtrPoints_%1.stl", 2 * m_stlLineWidth);

        writeOutput(QString(basename).append(QString(" - control polygons stored. Time: %1 ms").arg(timeStep->restart())));
    }

    if (doKnotLines) {
        const QVector<double> *knots = m_bSplineProperties->getKnotsP();
        const int lastKnotIndex = m_bSplineProperties->getLastInnerKnotIndex();
        QVector<double> knotlines;
        for (int i = m_bSplineProperties->getFirstInnerKnotIndex(); i <= lastKnotIndex; i++)  //only consider inner knot vectors
            knotlines << knots->at(i);

        m_stlWriter->writeBSplineIsolinesToStl(m_bSplinePatchNetwork->getBSplineSurfaces(), knotlines, knotlines, QString(basename).append("_knotlines_LW%1.stl").arg(m_stlLineWidth), m_parameterInterval, numberOfEvalPoints, m_stlLineWidth);
        writeOutput(QString(basename).append(QString(" - knot lines stored. Time %1 ms").arg(timeStep->restart())));
    }

    if (doDataError) {
        const double dataErrorEdgeWidth = 0;
        const bool doErrorOnBoundaries = true;

        if (doDataErrorForIndividualSurfaces)
            m_stlWriter->writeBSplineDataErrorToStl(m_bSplinePatchNetwork, QString(basename).append("_dataError.stl"), QString("debug/").append(basename).append("_dataError_%1.stl"), doErrorOnBoundaries, dataErrorEdgeWidth);
        else
            m_stlWriter->writeBSplineDataErrorToStl(m_bSplinePatchNetwork, QString(basename).append("_dataError.stl"), QString(), doErrorOnBoundaries, dataErrorEdgeWidth);

        writeOutput(QString(basename).append(QString(" - data error stored. Time: %1 ms").arg(timeStep->restart())));
    }
}

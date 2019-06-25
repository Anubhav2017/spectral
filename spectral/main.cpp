#include<QApplication>
#include <QCoreApplication>
#include <QString>
#include <QStringList>
#include <QInputDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include<laplacesolver.h>
#include <iostream>
#include <string>
#include <cstring>
#include "mainwindow.h"
#include "controller.h"
#include "stlwriter.h"
#include <QDebug>

#include"spectralsolver.cpp"

using namespace std;

int main(int argc, char *argv[]) {


    //QApplication a(argc, argv); //Needed for input widgets (like textbox)
    qDebug()<<"Hi";


    //get filename for the input stl
//    QString filename = QFileDialog::getOpenFileName(0, QString("Load triangulation"), QString("../.."), "*.stl");
//    if (filename.isEmpty())
//        return -1;

    //QString filename="/u/a/agarwala/Desktop/AdaptiveSurfaceReconstruction/project_files/stl_files/sphere-d050-sf1-ml_smoothed.stl";

    QString filename="/u/a/agarwala/Desktop/sphere-d050-sf1-ml_smoothed_downsampled.stl";

    //load stl file
    qDebug()<<"stl loaded";
    Mesh *m = Mesh::loadFromBinaryStl(filename);
    qDebug()<<"Mesh loaded";
    int ef=10;

    const int numberOfMeshFaces = m->getNumberOfFaces();
    qDebug() << "mesh loaded. Number of mesh faces" << m->getNumberOfFaces();

    QVector<double> v= computeSpectralEigenvector(m,ef);

    //ls.writeMeshWithVector(v,QString("morsefileTangle%1+_s").arg(ef));
    saveEigenvectorToTxt(m,v,QString("morsefiledownsample%1+_s").arg(ef));




//    QVector<QColor> colorMap=ls.generateColorMap2(vdash);
//    qDebug()<< "Color map generated";

//    ls.writeMeshWithVertexColors(colorMap,QString("Tanglevector%1s.ply").arg(ef));
//    qDebug()<< "color Mesh generated";

    QVector<int> vertexmap = loadMorseSmaleDecomposition("/u/a/agarwala/Desktop/spectral/quadmap10.txt");
 //   CellMesh* cm=ls.generateCellMesh("/u/a/agarwala/Desktop/spectral/quadmap10.txt");
    CellMesh* cm=createCellMeshFromMorseSmale(m,vertexmap);



    //test face and vertex valencies
    qDebug() << "Vertex Valencies" << cm->getValencyHistogramCellVertices() << "total:" << cm->getNumberOfVertices();
         qDebug() << "Face Valencies" << cm->getValencyHistogramCellFaces() << "total:" << cm->getNumberOfFaces();

    //Output of cellmesh into .stl
    QString outputFolder("../");
    QString basename("cellmesh_5");
    double stlLineWidth = 0.15;

    StlWriter writer(outputFolder, false);
    cm->cleanUpMeshIndices();   //Neccessary to avoid NULLs in vertex, edge, and face lists

    const int numberOfCellEdges = cm->getNumberOfEdges();
    QVector<Edge *> surfaceCellMeshEdges;
    for (int i = 0; i < numberOfCellEdges; i++) {
        if (cm->getEdge(i)) {
            QVector<PolylineVertex *> *polyVertices = cm->getEdge(i)->getPolylineVertices();
            for (int j = 1; j < polyVertices->size(); j++) {
                surfaceCellMeshEdges.append(new Edge(polyVertices->at(j-1), polyVertices->at(j), -1));
            }
        }
    }
    writer.writeEdgesToStl(&surfaceCellMeshEdges, QString(basename).append(QString("_borders_LW%1.stl").arg(stlLineWidth)), stlLineWidth);
    foreach (Edge *e, surfaceCellMeshEdges)
        delete e;
    qDebug() << "All done!";

    return 0;
}

//int main2(int argc, char *argv[])
//{
//    QCoreApplication a(argc, argv); //Needed for input widgets (like textbox)

//    const QString optionsIdentifierDefaultStlFolder("defaultStlFolder");
//    const QString optionsIdentifierSaveEnergyMatrix("saveEnergyMatrices");
//    const QString optionsIdentifierLoadEnergyMatrix("loadEnergyMatrices");
//    const QString optionsIdentifierEnergyFolder("energyMatrixFolder");
//    const QString optionsIdentifierLastInput("lastInput");

//    //if (argc = 1)
//    //  open GUI with most recent paramters
//    //if (argv containts -GUI)
//    //  open GUI with input paramters

//    //set up main classes
//    Controller algorithm;

//    //read options.ini
//    QVector<QString> optionIdentifiers;
//    QVector<QString> optionValues;
//    QFile optionsFile("./options.ini");
//    if (optionsFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
//        QTextStream optionsInStream(&optionsFile);
//        QString line;
//        while (!optionsInStream.atEnd()) {
//            line = optionsInStream.readLine();
//            QStringList lineSplit = line.split("=");
//            if (lineSplit.size() == 2) {
//                optionIdentifiers << lineSplit[0];
//                optionValues << lineSplit[1];
//            }
//        }
//    }
//    optionsFile.close();

//    //process options (and set up values if needed)
//    const int indexDefaultStlFolder = optionIdentifiers.indexOf(optionsIdentifierDefaultStlFolder);
//    QString stlFolder("../out/");
//    if (indexDefaultStlFolder >= 0)
//        stlFolder = optionValues[indexDefaultStlFolder];
//    else {
//        stlFolder = QFileDialog::getExistingDirectory(0, QString("Set Default .stl Output Folder"));
//        if (stlFolder.isEmpty()) {
//            std::cerr << "Error: No .stl Folder Selected!";
//            return -1;
//        }

//        optionIdentifiers << optionsIdentifierDefaultStlFolder;
//        optionValues << stlFolder;
//    }

//    const int indexSaveEnergyMatrix = optionIdentifiers.indexOf(optionsIdentifierSaveEnergyMatrix);
//    bool saveEnergyMatrix = false;
//    if (indexSaveEnergyMatrix >= 0)
//        saveEnergyMatrix = (optionValues[indexSaveEnergyMatrix].toInt() > 0);
//    else {
//        saveEnergyMatrix = (QMessageBox::question(0, "[Setup]", "Store Energy Matrices by Default?", QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes);
//        optionIdentifiers << optionsIdentifierSaveEnergyMatrix;
//        if (saveEnergyMatrix)
//            optionValues << "1";
//        else
//            optionValues << "0";

//        algorithm.setEnergyMatrixStorageOnHarddrive(saveEnergyMatrix);
//    }

//    const int indexLoadEnergyMatrix = optionIdentifiers.indexOf(optionsIdentifierLoadEnergyMatrix);
//    bool loadEnergyMatrix = false;
//    if (indexLoadEnergyMatrix >= 0)
//        loadEnergyMatrix = (optionValues[indexLoadEnergyMatrix].toInt() > 0);
//    else {
//        loadEnergyMatrix = (QMessageBox::question(0, "[Setup]", "Load Energy Matrices by Default?", QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes);
//        optionIdentifiers << optionsIdentifierLoadEnergyMatrix;
//        if (loadEnergyMatrix)
//            optionValues << "1";
//        else
//            optionValues << "0";

//        algorithm.setEnergyMatrixLoadingFromHarddrive(loadEnergyMatrix);
//    }

//    if (saveEnergyMatrix || loadEnergyMatrix) {
//        const int indexEnergyFolder = optionIdentifiers.indexOf(optionsIdentifierEnergyFolder);
//        if (indexEnergyFolder >= 0)
//            algorithm.setEnergyFolder(optionValues[indexEnergyFolder]);
//        else {
//            QString energyFolder = QFileDialog::getExistingDirectory(0, QString("Select Energy Matrix Storage Location"));
//            if (energyFolder.isEmpty()) {
//                std::cerr << "Error: No Folder for Energy Storage Selected!";
//                return -1;
//            }

//            optionIdentifiers << optionsIdentifierEnergyFolder;
//            optionValues << energyFolder;
//            algorithm.setEnergyFolder(energyFolder);
//        }
//    }

//    QString lastInput;
//    const int indexLastInput = optionIdentifiers.indexOf(optionsIdentifierLastInput);
//    if (indexLastInput >= 0)
//        lastInput = optionValues[indexLastInput];

//    //prepare input arguments
//    int numArgs = 0;
//    QStringList arguments;
//    if (argc == 1) {
//        //if no command line arguments are given, use text dialog to ask for arguments
//        bool okPressed;
//        QString text = QInputDialog::getText(0, QString("Input Arguments"), QString("Enter parameters: "), QLineEdit::Normal, lastInput, &okPressed);
//        if (!okPressed)
//            exit(0);
//        if (!text.isEmpty()) {
//            arguments = text.split(" ");
//            numArgs = arguments.size();
//        }

//        //save text input
//        if (indexLastInput >= 0)
//            optionValues[indexLastInput] = text;
//        else {
//            optionIdentifiers << optionsIdentifierLastInput;
//            optionValues << text;
//        }
//    } else {
//        numArgs = argc - 1;
//        for (int i = 1; i < argc; i++)
//            arguments.append(QString(argv[i]));
//    }

//    //Rewrite options file
//    if (optionsFile.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
//        QTextStream optionsOutStream(&optionsFile);
//        const int numLines = optionIdentifiers.size();
//        for (int i = 0; i < numLines; i++)
//            optionsOutStream << optionIdentifiers[i] << "=" << optionValues[i] << "\n";
//    }
//    optionsFile.close();

//    QString inputArguments = QString(argv[0]).split("\\").last().split("/").last();
//    for (int i = 0; i < numArgs; i++)
//        inputArguments.append(QString(" ") + arguments[i]);

//    //parse input arguments
//    int i = 0;
//    while (i < numArgs) {
//        QString arg = arguments[i];
//        if (arg == "-help" || arg == "-h" || arg == "-?") {
//            std::cout << "Usage: " << std::endl
//                      << " -help | -h | -?   show help" << std::endl
//                      << "General Parameters:" << std::endl
//                      << " -s [int]          set start step in algorithm" << std::endl
//                      << " -e [int]          set end step in algorithm" << std::endl
//                      << " -dm [string]      set folder for results" << std::endl
//                      << " -stl [string]     set folder for stl export" << std::endl
//                      << " -stlRes [int]     set number of sample points for converting B-splines to .stl" << std::endl
//                      << " -stlLW [double]   set line width for curves and polylines in .stl files" << std::endl
//                      << " -d [int]          set debug mode" << std::endl
//                      << " -noConsole        deactivate console output" << std::endl
//                      << " -noLog            deactivate log output" << std::endl
//                      << " -clearlog         clear the log file before starting" << std::endl
//                      << " -t [int]          set number of threads (disable automatic thread count)" << std::endl
//                      << "Step 0 - obtain mesh" << std::endl
//                      << " -scaleMesh [int]  scale mesh to fit into (0) unit cube (1) sphere with d=1" << std::endl
//                      << " -keepProportions  when scaling, keep the proportions of the mesh" << std::endl
//                      << " -sf [int]         scaling factor: edge length of target cube" << std::endl
//                      << "Step 1 - smooth mesh" << std::endl
//                      << " -smoothIt [int]   set (max.) number of iterations for mesh smoothing" << std::endl
//                      << "Step 2 - compute quadmesh" << std::endl
//                      << " -nic [int]        set minimal number of cells in initial (Voronoi) tessellation" << std::endl
//                      << " -nqc [int]        set minimal number of cells in quadmesh" << std::endl
//                      << " -qma [d,v,c]      set quadmesh algorithm to [d] DelaunayMatching, [v] VoronoiSubdivision, [c] CuttingPlanes (manual)" << std::endl
//                      << " -instantTopo      do instant topology check in Voronoi algorithm" << std::endl
//                      << " -skipRef          skip the topology improving refinement phase" << std::endl
//                      << " -qmsd [int]       set number of (additional) central subdivisions for quadmesh" << std::endl
//                      << " -param [double] [double] set interval for paramterization" << std::endl
//                      << "Step 3 - optimize quadmesh" << std::endl
//                      << " -qmo [int]        set number of vertex optimizations" << std::endl
//                      << "Step 4 - do B-spline fitting" << std::endl
//                      << " -bsa [cg0, cg0i, gg0, sgg1, lg1] set BSpline algorithm" << std::endl
//                      << " -invertQM         invert orientation of each cell in the quadmesh" << std::endl
//                      << " -newParam         recompute the parameterization of the quadmesh (harmonic map)" << std::endl
//                      << " -initGuess        use existing results as initial guess" << std::endl
//                      << " -bnc [int]        set number of B-spline control points per row" << std::endl
//                      << " -bo [int]         set spline order" << std::endl
//                      << " -bwe [double]     set weight for energy term in B-spline fitting (LS=1-Energy)" << std::endl
//                      << " -bvg1r [double]   set parameter range of vertex G1 fitting" << std::endl
//                      << " -bvg1n [int]      set number of G1 points for local vertex fitting" << std::endl
//                      << " -bvg1Opt          do vertex fitting with G1 optimization (no G1 constraints)" << std::endl
//                      << " -bvg1w [double]   set weight for G1 constraints in vertex optimization" << std::endl
//                      << " -beg1n [int]      set number of G1 points for local edge fitting" << std::endl
//                      << " -beg1Opt          do edge fitting with G1 optimization (no G1 constraints)" << std::endl
//                      << " -beg1w [double]   set weight for G1 constraints in edge optimization" << std::endl
//                      << " -detFit           activate detailed output for spline fitting" << std::endl
//                      << "Step 5 - optimize paramter values" << std::endl
//                      << " -nopo             do not optimize parameters after fitting" << std::endl
//                      << " -po [int]         set number of parameter optimizations after the fitting step" << std::endl
//                      << " -pnlo             use nonlinear optimization to find optimal parameters" << std::endl
//                      << "Step 6 - evaluate B-splines" << std::endl
//                      << " -scaleEval [int]  use scaling for evaluation (parameters are the same as for mesh scaling above)" << std::endl
//                      << " -noMatlab         do not write matlab output files" << std::endl
//                      << " -do3Dimage        do 3D image output" << std::endl
//                      << " -3DSize [int, int, int]  set size of resulting 3D image (otherwise, it will use the bounding box)" << std::endl
//                      << "Step 7 - model based surface inspection" << std::endl
//                      << " -mbsiMet [int]      set metric type (0 - TP, 1 - Curvature, 2 - Area, 3 - Normal deviation, 4 - Inspection Angle, >= 5 - Highlight)" << std::endl
//                      << " -mbsiTh [double] [a, f] threshold for curvature based subdivision, a - meant as absolute value, f - factor of average metric value" << std::endl
//                      << " -mbsiN [int]        instead of threshold based subdivision, produce N viewpoint candidates" << std::endl
//                      << " -mbsiRef            do reference algorithms as well (every vertex and random vertex)" << std::endl
//                      << " -mbsiCDist [double] set distance of camera to pivot point" << std::endl
//                      << " -mbsiPSize [double] set size of pivot point (for stl output)" << std::endl
//                      << " -mbsiCSize [double] set size of camera point (for stl output)" << std::endl
//                      << " -mbsiVSize [double] set size of direction vector (for stl output)" << std::endl
//                      << " -mbsiVL [double]    set length of direction vector (for stl output, default = cam distance)" << std::endl
//                      << "Parameters for numerics:" << std::endl
//                      << " -nloIt [int]      set number of iterations for nlopt solvers" << std::endl
//                      << " -nloFP [double]   set relative abort criterion for objective function" <<  std::endl
//                      << " -nloXP [double]   set relative abort criterion for objective variables" << std::endl
//                      << " -nloCP [double]   set precision for equality constraints" << std::endl;

//            return 0;
//        } else if (arg == "-s") {              //start point
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setStartPoint(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-e") {       //end point
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setEndPoint(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-dm") {      //dm folder
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            QString folder(arguments[i + 1]);
//            if (!folder.endsWith("/"))
//                folder.append("/");
//            algorithm.setDmFolder(folder);

//            i += 2;
//        } else if (arg == "-stl") {     //stl folder
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            stlFolder = arguments[i + 1];

//            i += 2;
//        } else if (arg == "-stlRes") {  //set number of sampling points for conversion of B-splines to stl
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setStlBSplineSamplePoints(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-stlLW") {   //line width for curves and polylines
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setStlLineWidth(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-d") {       //debuging mode
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setDebugingMode(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-noConsole") {
//            algorithm.setConsoleOutput(false);

//            i += 1;
//        } else if (arg == "-noLog") {
//            algorithm.setLogOutput(false);

//            i += 1;
//        } else if (arg == "-clearlog") {//clear log file
//            algorithm.clearLogFile();

//            i += 1;
//        } else if (arg == "-t") {       //number of threads
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfThreads(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-scaleMesh") {   //scale mesh into unit cube (unless scaling factor is specified)
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMeshScaling(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-keepProportions") { //keep mesh proportions when scaling
//            algorithm.setKeepingProportions(true);

//            i += 1;
//        } else if (arg == "-sf") {  //multiply scaled mesh by a scaling factor (default: 1)
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setScalingFactor(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-smoothIt") {//max number of smoothing iterations
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setmaxSmoothIterations(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-nic") {     //min number of initial cells in voronoi tessellation
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMinNumberOfInitialCells(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-nqc") {     //min number of cells in quadmesh
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMinNumberOfQuadCells(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-qma") {     //set quadmesh algorithm
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            if (arguments[i + 1] == "d")
//                algorithm.setQuadmeshAlgorithm(Controller::QMA_DelaunayMatching);
//            else if (arguments[i + 1] == "v")
//                algorithm.setQuadmeshAlgorithm(Controller::QMA_VoronoiSubdivision);
//            else if (arguments[i + 1] == "c")
//                algorithm.setQuadmeshAlgorithm(Controller::QMA_CuttingPlanes);
//            else {
//                std::cerr << "unknown algorithm identifier " << arguments[i + 1].toStdString() << std::endl;
//                return -1;
//            }

//            i += 2;
//        } else if (arg == "-instantTopo") { //activate instant topology check (for voronoi algorithms)
//            algorithm.setInstantTopologyCheck(true);
//            i += 1;
//        } else if (arg == "-skipRef") {
//            algorithm.setSkipRefinement(true);
//            i += 1;
//        } else if (arg == "-qmsd") {     //number of central subdivisions for quadmesh
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfCentralSubdivisions(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-param") {   //set parameter interval
//            if (i + 2 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setParameterInterval(arguments[i + 1].toDouble(), arguments[i + 2].toDouble());

//            i += 3;
//        } else if (arg == "-qmo") {     //number of vertex optimizations in quadmesh
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfVertexOptimizations(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-bsa") {
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            if (arguments[i + 1] == "cg0")
//                algorithm.setBSplineAlgorithm(Controller::BSA_CurveFittingG0);
//            else if (arguments[i + 1] == "cg0i")
//                algorithm.setBSplineAlgorithm(Controller::BSA_CurveFittingG0Iterative);
//            else if (arguments[i + 1] == "gg0")
//                algorithm.setBSplineAlgorithm(Controller::BSA_GlobalG0);
//            else if (arguments[i + 1] == "sgg1")
//                algorithm.setBSplineAlgorithm(Controller::BSA_SimplifiedGlobalG1);
//            else if (arguments[i + 1] == "lg1")
//                algorithm.setBSplineAlgorithm(Controller::BSA_LocalG1Multiphase);
//            else {
//                std::cerr << "unknown algorithm identifier " << arguments[i + 1].toStdString() << std::endl;
//                return -1;
//            }

//            i += 2;
//        } else if (arg == "-invertQM") {//Invert orientation of quadmesh
//            algorithm.setInvertCellmeshOrientation(true);
//            i += 1;
//        } else if (arg == "-newParam") {//recompute parameterization
//            algorithm.setRecomputeParameterization(true);
//            i += 1;
//        } else if (arg == "-initGuess") {//use existing B-spline network as initial guess
//            algorithm.setUsageOfInitialGuess(true);
//            i += 1;
//        } else if (arg == "-bnc") {     //number of b-spline control points per row
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfControlPoints(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-bo") {      //b-spline order
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setSplineOrder(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-bwe") {     //weight for energy term
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setWeightEnergyTerm(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-bvg1r") {   //vertex g1 range
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setG1VertexParameterRange(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-bvg1n") {   //vertex g1 points
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfG1VertexPoints(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-bvg1Opt") { //use g1 constraints as objective values for nlo (vertex)
//            algorithm.setVertexFittingG1Opt(true);
//            i += 1;
//        } else if (arg == "-bvg1w") {   //weight of g1 constraints in objective function (vertex)
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setVertexFitG1Weight(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-beg1n") {   //edge g1 points
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfG1EdgePoints(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-beg1Opt") { //use g1 constraints as objective values for nlo (edge)
//            algorithm.setEdgeFittingG1Opt(true);
//            i += 1;
//        } else if (arg == "-beg1w") {   //weight of g1 constraints in objective function (edge)
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setEdgeFitG1Weight(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-detFit") {  //activate detailed output for spline fitting
//            algorithm.setDetailedFittingOutput(true);
//            i += 1;
//        } else if (arg == "-nopo") {    //do not optimize parameters
//            algorithm.setNumberOfParameterOptimizations(0);
//            i += 1;
//        } else if (arg == "-po") {      //number of parameter optimizations
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNumberOfParameterOptimizations(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-pnlo") {    //use nonlinear optimization for parameter optimization
//            algorithm.setNumberOfParameterOptimizations(1);
//            algorithm.setUseOfNLOforParamOpt(true);

//            i += 1;
//        } else if (arg == "-scaleEval") {   //scale mesh and bsplines for evaluation only
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setScalingForEvaluation(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-noMatlab") {    //do not write matlab files
//            algorithm.setMatlabOutput(false);

//            i += 1;
//        } else if (arg == "-do3Dimage") {   //do 3D image output
//            algorithm.set3DImageOutput(true);

//            i += 1;
//        } else if (arg == "-3DSize") {  //set size of 3D output image
//            if (i + 3 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.set3DImageOutput(true);   //currently there is no other reason to set 3D image sizes
//            algorithm.set3DImageSizes(arguments[i + 1].toInt(), arguments[i + 2].toInt(), arguments[i + 3].toInt());

//            i += 4;
//        } else if (arg == "-mbsiMet") { //metric type
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiMetricType(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-mbsiTh") {  //threshold for subdivision
//            if (i + 2 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            const bool avgFactor = (arguments[i + 2] == "f");

//            if (arguments[i + 2] != "f" && arguments[i + 2] != "a") {
//                std::cerr << "Unknown identifier to set mbsi threshold: " << arguments[i + 2].toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiThreshold(arguments[i + 1].toDouble(), avgFactor);

//            i += 3;
//        } else if (arg == "-mbsiN") {   //set target number of VPC
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiTargetNumberOfVPC(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-mbsiRef") { //Use reference algorithms, too
//            algorithm.setMbsiDoReferenceAlgorithms(true);

//            i += 1;
//        } else if (arg == "-mbsiCDist") {   //camera distance
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiCamDistance(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-mbsiPSize") {   //for stl output: size of pivot points
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiPivotPointSize(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-mbsiCSize") {   //for stl output: size of camera positions
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiCamPointSize(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-mbsiVSize") {   //for stl output: size of vectors between pivot points and cameras
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiCamDirectionSize(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-mbsiVL") {      //for stl output: length of normal vector
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setMbsiCamDirectionLength(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-nloIt") {
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNloptMaxIterations(arguments[i + 1].toInt());

//            i += 2;
//        } else if (arg == "-nloFP") {
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNloptFPrec(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-nloXP") {
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNloptXPrec(arguments[i + 1].toDouble());

//            i += 2;
//        } else if (arg == "-nloCP") {
//            if (i + 1 >= numArgs) {
//                std::cerr << "Missing command line arguments after " << arg.toStdString() << std::endl;
//                return -1;
//            }

//            algorithm.setNloptConstrPrec(arguments[i + 1].toDouble());

//            i += 2;
//        } else {
//            std::cerr << "Invalid command line argument: " << arg.toStdString() << std::endl;
//            return -1;
//        }
//    }

//    if (!stlFolder.endsWith("/"))
//        stlFolder.append("/");

//    StlWriter debugWriter(stlFolder, false, false); //TODO rethink this concept (especially when the spline generator gets its own thread)
//    algorithm.setStlWriter(&debugWriter);

//    bool showGUI = false;
//    if (showGUI) {
//        MainWindow w;
//        w.show();
//    }

//    algorithm.writeToLogOnly(QString());    //Adds a newline to the log
//    algorithm.writeOutput(inputArguments);

//    algorithm.start();

//    if (showGUI)
//        return a.exec();
//    else
//        return 0;
//}

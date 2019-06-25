#-------------------------------------------------
#
# Project created by QtCreator 2017-10-16T14:33:39
#
#-------------------------------------------------
QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = spectral
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app



SOURCES += main.cpp\
    mainwindow.cpp \
    Mesh/edge.cpp \
    Mesh/vertex.cpp \
    Mesh/face.cpp \
    Mesh/mesh.cpp \
    CellMesh/cellvertex.cpp \
    CellMesh/celledge.cpp \
    CellMesh/cellface.cpp \
    CellMesh/cellmesh.cpp \
    CellMesh/polylinevertex.cpp \
    CellMesh/parameterization.cpp \
    SurfaceTessellation/surfacevoronoigenerator.cpp \
    SurfaceTessellation/quadmeshgenerator.cpp \
    SurfaceTessellation/blossomgraphmatching.cpp \
    BSpline/bsplinesurface.cpp \
    BSpline/bsplinecurve.cpp \
    BSpline/bsplinepatchnetwork.cpp \
    BSpline/bsplinetriangulation.cpp \
    MBSI/viewpointcandidategenerator.cpp \
    MBSI/subdivisionthinplateenergymetric.cpp \
    MBSI/subdivisionexactcurvaturemetric.cpp \
    MBSI/subdivisionsurfaceareametric.cpp \
    MBSI/subdivisionmanualhighlightmetric.cpp \
    MBSI/subdivisionnormaldeviationmetric.cpp \
    MBSI/colormapgenerator.cpp \
    stlwriter.cpp \
    splinefitter.cpp \
    parameteroptimizer.cpp \
    voxelizer.cpp \
    splineevaluator.cpp \
    external/DataIO.cpp \
    external/IASS.cpp \
    controller.cpp \
    energymatrixgenerator.cpp \
    plywriter.cpp \
    BSpline/bsplineproperties.cpp \
    MBSI/surfacesubdivider.cpp \
    MBSI/subsurfacetree.cpp \
    MBSI/subdivisioninspectionanglemetric.cpp \
    laplacesolver.cpp \
    spectralsolver.cpp

HEADERS  += mainwindow.h\
    Mesh/edge.h \
    Mesh/vertex.h \
    Mesh/face.h \
    Mesh/mesh.h \
    CellMesh/cellvertex.h \
    CellMesh/celledge.h \
    CellMesh/cellface.h \
    CellMesh/cellmesh.h \
    CellMesh/polylinevertex.h \
    CellMesh/parameterization.h \
    SurfaceTessellation/priorityqueue.h \
    SurfaceTessellation/surfacevoronoigenerator.h \
    SurfaceTessellation/quadmeshgenerator.h \
    SurfaceTessellation/blossomgraphmatching.h \
    BSpline/bsplinesurface.h \
    BSpline/bsplinecurve.h \
    BSpline/bsplinepatchnetwork.h \
    BSpline/bsplinetriangulation.h \
    MBSI/viewpointcandidategenerator.h \
    MBSI/subdivisionabstractmetric.h \
    MBSI/subdivisionthinplateenergymetric.h \
    MBSI/subdivisionexactcurvaturemetric.h \
    MBSI/subdivisionsurfaceareametric.h \
    MBSI/subdivisionmanualhighlightmetric.h \
    MBSI/subdivisionnormaldeviationmetric.h \
    MBSI/colormapgenerator.h \
    stlwriter.h \
    splinefitter.h \
    parameteroptimizer.h \
    voxelizer.h \
    splineevaluator.h \
    external/DataIO.h \
    external/Miniball.hpp \
    external/IASS.h \
    controller.h \
    energymatrixgenerator.h \
    plywriter.h \
    BSpline/bsplineproperties.h \
    MBSI/surfacesubdivider.h \
    MBSI/subsurfacetree.h \
    MBSI/subdivisioninspectionanglemetric.h \
    quadmeshgeneration.h \
    laplacesolver.h

FORMS    += mainwindow.ui

# activate open mp
QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

# Eigen v.3.3.2 [2017_02_05]
# http://eigen.tuxfamily.org/index.php?title=Main_Page
INCLUDEPATH += "../external/Eigen"
INCLUDEPATH += "../external"


INCLUDEPATH += "../external/igl"

# GNU Scientific Library (GSL) [2017_02_05]
# -> Linux (v.2.3): https://www.gnu.org/software/gsl/
# -> Windows (v.1.8): http://gnuwin32.sourceforge.net/packages/gsl.htm

INCLUDEPATH += "../external/gsl/include"
LIBS += "../external/gsl/lib/libgsl.a"
LIBS += "../external/gsl/lib/libgslcblas.a"

# NLOpt 2.4.2 [2017_05_28]
# http://ab-initio.mit.edu/wiki/index.php/NLopt
INCLUDEPATH += "../external/nlopt-2.4.2/include"
win32 {
   LIBS += -L"../external/nlopt-2.4.2/lib/" -llibnlopt-0
}
unix {
   LIBS += "../external/nlopt-2.4.2/lib/libnlopt.a"
}

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/release/ -lCGAL
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/debug/ -lCGAL
else:unix: LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/ -lCGAL

INCLUDEPATH += $$PWD/../../CGAL-4.13.1/build/include
INCLUDEPATH += $$PWD/../../CGAL-4.13.1/include
DEPENDPATH += $$PWD/../../CGAL-4.13.1/build/include




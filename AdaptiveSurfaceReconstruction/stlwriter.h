#ifndef STLWRITER_H
#define STLWRITER_H

#include "Mesh/mesh.h"
#include "CellMesh/cellface.h"
#include "CellMesh/celledge.h"
#include "BSpline/bsplinesurface.h"
#include "BSpline/bsplinepatchnetwork.h"
#include "MBSI/subsurfacetree.h"

#include <QString>
#include <QFile>
#include <QVector3D>
#include <QTextStream>

class StlWriter
{
public:
    StlWriter(QString defaultFolder = QString(), bool cleanDirectoryBeforeWriting = true, bool askOverwrite = false, bool consoleOutput = false);
    void copyAllFilesTo(QString destinationFolder);

    void writeMeshFacesToStl(Mesh *m, QString filename);
    void writeFacesToStlBinary(QVector<Face *> *faces, QString filename);
    void writeFacesToStlASCII(QVector<Face *> *faces, QString filename);
    void writeEdgesToStl(QVector<Edge *> *edges, QString filename, double width=0);
    void writeBSplineCurvesToStl(QVector<BSplineCurve *> curves, QString filename, double linewidth = 0.2, int numberOfEvalPointsPerCurve = 51);
    void writeBSplineSurfacesToStl(QVector<BSplineSurface *> surfaces, QString filename, int numberOfEvalPointsPerRow = 51);
    void writeBSplineSubdivisionsToStl(QVector<BSplineSurface *> surfaces, QVector<SubSurfaceTree *> subdivisions, QString filename, double linewidth = 0, int numberOfEvalPoints = 51);
    void writeCellFaceContentToStl(CellFace *cf, QString filename, double vertexSize = 0.2, double edgeWidth = 0);
    void writeCellBorderContentToStl(CellFace *cf, Mesh *originalMesh, QString filename, double vertexSize = 0.2, double edgeWidth = 0);
    void writeControlPointsOfIndividualPhases(BSplinePatchNetwork *patchNetwork, QString filenameTemplate, double vertexSize = 0.2);
    void writeBSplineControlPointsToStl(BSplinePatchNetwork *patchNetwork, QString filenameEdges, QString filenameVertices, double vertexSize = 0.2, double linewidth = 0.2);
    void writeBSplineIsolinesToStl(QVector<BSplineSurface *> surfaces, QVector<double> isolinesU, QVector<double> isolinesV, QString filename, QPair<double, double> parameterInterval, int numberOfEvalPoints = 51, double linewidth = 0.2);
    void writeBSplineDataErrorToStl(BSplinePatchNetwork *patchNetwork, QString filename, QString indSurfaceNameTemplate, bool errorOnBoundaries = false, double edgeWidth = 0);

    void writePoints(QVector<QVector3D> *points, double size, QString filename);

    void setAskOverwrite(bool newValue);
    bool getAskOverwrite();

    void setConsoleOutput(bool newValue);
    bool getConsoleOutput();

    void setDefaultFolder(QString defaultFolder);
    QString getDefaultFolder();
private:
    QString m_defaultFolder;
    bool m_askOverwrite;
    bool m_consoleOutput;

    void writeLineSegment(const QVector3D &start, const QVector3D &end, QTextStream *out, const double width);
    QFile *processFilename(QString filename);
    bool removeRecursively(const QString &dirName);
    bool copyRecursively(const QString &srcFilePath, const QString &tgtFilePath);
};

#endif // STLWRITER_H

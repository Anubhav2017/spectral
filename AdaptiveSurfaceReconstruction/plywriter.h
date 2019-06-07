#ifndef PLYWRITER_H
#define PLYWRITER_H

#include "BSpline/bsplinepatchnetwork.h"
#include "BSpline/bsplinetriangulation.h"
#include "Mesh/mesh.h"
#include <QVector>
#include <QColor>

class PlyWriter
{
public:
    static void writeBSplineTriangulationWithVertexColors(BSplinePatchNetwork *bSplines, BSplineTriangulation *triangulation, const QVector<QColor> &colormap, const QString filename);
    static void writeMeshWithVertexColors(Mesh *m, const QVector<QColor> &colormap, const QString filename);
};

#endif // PLYWRITER_H

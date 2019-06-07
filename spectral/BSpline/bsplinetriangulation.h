#ifndef BSPLINETRIANGULATION_H
#define BSPLINETRIANGULATION_H

#include <QVector>
#include <QVector2D>

class BSplineTriangulation
{
public:
    static BSplineTriangulation *createRegularSampledTriangulation(const int numSurfaces, const int numSamplePoints, const double u0, const double u1, const double v0, const double v1, const bool quadsInsteadOfTris);

    int getNumberOfSamplePointsPerSurface();
    int getNumberOfSamplePointsPerRow();
    int getNumberOfSurfaces();
    int getNumberOfVertices();
    int getNumberOfFaces();
    int getSurfaceIdOfVertex(int id);
    int getVertexOffsetForSurface(int surfaceId);
    QVector2D getParameterValueForVertex(int id);
    QVector<int> getFace(int id);
private:
    BSplineTriangulation(const int numSurfaces, const int numSamplePoints, const QVector<QVector2D> &paramValues, const QVector<QVector<int> > &faces);
    const int m_numberOfSurfaces;
    const int m_numberOfSamplePoints;

    QVector<QVector2D> m_vertexParameterValues;
    QVector<QVector<int> > m_faces;
};

#endif // BSPLINETRIANGULATION_H

#include "bsplinetriangulation.h"

BSplineTriangulation *BSplineTriangulation::createRegularSampledTriangulation(const int numSurfaces, const int numSamplePointsPerRow, const double u0, const double u1, const double v0, const double v1, const bool quadsInsteadOfTris)
{
    //Sample point indexing
    //
    // ...
    // ...
    // n+1 n+2 n+3
    // 1   2   3    ...  n

    const int numPointsPerSurf = numSamplePointsPerRow * numSamplePointsPerRow;
    const double uDiff = u1 - u0;
    const double vDiff = v1 - v0;

    QVector<QVector<int> > faces;
    QVector<QVector2D> paramValues(numPointsPerSurf * numSurfaces);

    for (int iSurf = 0; iSurf < numSurfaces; iSurf++) {
        const int offset = iSurf * numPointsPerSurf;

        for (int iu = 0; iu < numSamplePointsPerRow; iu++) {
            for (int iv = 0; iv < numSamplePointsPerRow; iv++) {
                const int indexLowLeft = offset + (iv * numSamplePointsPerRow) + iu;
                const double u = u0 + uDiff * ((double) iu / (numSamplePointsPerRow - 1));
                const double v = v0 + vDiff * ((double) iv / (numSamplePointsPerRow - 1));
                paramValues[indexLowLeft] = QVector2D(u, v);

                if (iu < numSamplePointsPerRow - 1 && iv < numSamplePointsPerRow - 1) {
                    const int indexLowRight = indexLowLeft + 1;
                    const int indexUpLeft = indexLowLeft + numSamplePointsPerRow;
                    const int indexUpRight = indexUpLeft + 1;
                    if (quadsInsteadOfTris) {
                        QVector<int> quad;
                        quad << indexLowLeft << indexLowRight << indexUpRight << indexUpLeft;
                        faces << quad;
                    } else {
                        QVector<int> tri0;
                        tri0 << indexLowLeft << indexLowRight << indexUpRight;
                        QVector<int> tri1;
                        tri1 << indexLowLeft << indexUpRight << indexUpLeft;
                        faces << tri0 << tri1;
                    }
                }
            }
        }
    }

    return new BSplineTriangulation(numSurfaces, numSamplePointsPerRow, paramValues, faces);
}

int BSplineTriangulation::getNumberOfSamplePointsPerSurface()
{
    return m_numberOfSamplePoints * m_numberOfSamplePoints;
}

int BSplineTriangulation::getNumberOfSamplePointsPerRow()
{
    return m_numberOfSamplePoints;
}

int BSplineTriangulation::getNumberOfSurfaces()
{
    return m_numberOfSurfaces;
}

int BSplineTriangulation::getNumberOfVertices()
{
    return m_vertexParameterValues.size();
}

int BSplineTriangulation::getNumberOfFaces()
{
    return m_faces.size();
}

int BSplineTriangulation::getSurfaceIdOfVertex(int id)
{
    return id / (m_numberOfSamplePoints * m_numberOfSamplePoints);
}

int BSplineTriangulation::getVertexOffsetForSurface(int surfaceId)
{
    return m_numberOfSamplePoints * m_numberOfSamplePoints * surfaceId;
}

QVector2D BSplineTriangulation::getParameterValueForVertex(int id)
{
    return m_vertexParameterValues[id];
}

QVector<int> BSplineTriangulation::getFace(int id)
{
    return m_faces[id];
}

BSplineTriangulation::BSplineTriangulation(const int numSurfaces, const int numSamplePoints, const QVector<QVector2D> &paramValues, const QVector<QVector<int> > &faces):
    m_numberOfSurfaces(numSurfaces), m_numberOfSamplePoints(numSamplePoints),
    m_vertexParameterValues(paramValues), m_faces(faces)
{

}

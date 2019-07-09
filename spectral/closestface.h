#ifndef CLOSESTFACE_H
#define CLOSESTFACE_H
#include<QVector>
#include<QVector3D>
#include<Eigen>
#include<Mesh/mesh.h>
#include<igl/point_mesh_squared_distance.h>
#include<QDebug>
#include<stlwriter.h>
#include <igl/readOBJ.h>
class closestface
{
public:

    static QVector<int> loadFaceMap(Mesh* m1, QVector<QVector3D> pointsvector);
    static QVector<int> loadFaceMapTriMesh(Mesh* m1, QVector<QVector3D> pointsvector);
    static QVector<int> loadFaceMapQuadMesh(Mesh* m1, QVector<QVector3D> pointsvector);

    static QVector<QVector<double> > computeBarycentreCoordinates(Mesh* m, QVector<QVector3D> pointsvector);
    static QVector<QVector<double> > computeBarycentreCoordinatesTriMesh(Mesh* m, QVector<QVector3D> pointsvector);
    static QVector<QVector<double> > computeBarycentreCoordinatesQuadMesh(Mesh* m, QVector<QVector3D> pointsvector);

    static QVector<QVector3D> loadProjections(Mesh *m, QVector<QVector3D> pointsvector);
    static QVector<QVector3D> loadProjectionsTriMesh(Mesh *m, QVector<QVector3D> pointsvector);
    static QVector<QVector3D> loadProjectionsQuadMesh(Mesh *m, QVector<QVector3D> pointsvector);

    static QVector<QVector2D > parametrizationQuadMesh(Mesh* mquad, QVector<QVector3D> pointsvector);
};

#endif // CLOSESTFACE_H

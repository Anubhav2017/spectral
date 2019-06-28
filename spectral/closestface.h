#ifndef CLOSESTFACE_H
#define CLOSESTFACE_H
#include<QVector>
#include<QVector3D>
#include<Eigen>
#include<Mesh/mesh.h>
#include<igl/point_mesh_squared_distance.h>
#include<QDebug>
class closestface
{
public:

    static QVector<int> loadFaceMap(Mesh* m1, Eigen::MatrixXd points);
    static QVector<QVector<double> > computeBarycentreCoordinates(Mesh* m, Eigen::MatrixXd points);
    static QVector<QVector3D> loadProjections(Mesh *m, Eigen::MatrixXd points);
};

#endif // CLOSESTFACE_H

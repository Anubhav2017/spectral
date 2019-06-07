#ifndef BSPLINEPATCHNETWORK_H
#define BSPLINEPATCHNETWORK_H

#include "CellMesh/cellmesh.h"
#include "bsplinesurface.h"
#include "bsplinecurve.h"
#include "bsplineproperties.h"

class BSplinePatchNetwork
{
public:
    BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties);
    BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties, QVector<BSplineSurface *> &bSplineSurfaces);
    BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties, QVector<BSplineSurface *> &bSplineSurfaces, QVector<BSplineCurve *> &bSplineCurves);
    ~BSplinePatchNetwork();

    void getBoundingBox(int &minX, int &maxX, int &minY, int &maxY, int &minZ, int &maxZ);
    void translate(const QVector3D &translationVector);
    void scale(double factorX, double factorY, double factorZ);

    CellMesh *getCellMesh();
    int getNumberOfBSplineSurfaces();
    int getNumberOfBSplineCurves();
    BSplineProperties *getProperties();

    BSplineSurface *getBSplineSurface(int id);
    QVector<BSplineSurface *> getBSplineSurfaces(); //Return a copy of the vector (contains pointers to the same surfaces)
    BSplineCurve *getBSplineCurve(int id);
    QVector<BSplineCurve *> getBSplineCurves();     //Return a copy of the vector (contains pointers to the same curves)

    void setBSplines(QVector<BSplineSurface *> &bSplineSurfaces, QVector<BSplineCurve *> &bSplineCurves);
    void setBSplines(QVector<BSplineSurface *> &bSplineSurfaces);

    void addInfo(QString info);
    void saveToFile(QString filename, QString optionalHeader = QString());
    static BSplinePatchNetwork *loadFromFile(QString filename, CellMesh *cm);
private:
    void extractBorderCurves();
    CellMesh *m_cellMesh;

    BSplineProperties *m_properties;
    QVector<BSplineSurface *> m_bSplineSurfaces;
    QVector<BSplineCurve *> m_bSplineCurves;

    QString m_info;
};

#endif // BSPLINEPATCHNETWORK_H

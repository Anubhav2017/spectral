#include "bsplinepatchnetwork.h"

//Qt includes
#include <QVector>
#include <QVector3D>
#include <QFile>
#include <QTextStream>
#include <QStringList>

BSplinePatchNetwork::BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties) :
    m_cellMesh(cm), m_properties(properties)
{

}

BSplinePatchNetwork::BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties, QVector<BSplineSurface *> &bSplineSurfaces):
    m_cellMesh(cm), m_properties(properties), m_bSplineSurfaces(bSplineSurfaces)
{
    this->extractBorderCurves();
}

BSplinePatchNetwork::BSplinePatchNetwork(CellMesh *cm, BSplineProperties *properties, QVector<BSplineSurface *> &bSplineSurfaces, QVector<BSplineCurve *> &bSplineCurves):
    m_cellMesh(cm), m_properties(properties), m_bSplineSurfaces(bSplineSurfaces), m_bSplineCurves(bSplineCurves)
{

}

BSplinePatchNetwork::~BSplinePatchNetwork()
{
    const int numberOfSurfaces = m_bSplineSurfaces.size();
    const int numberOfCurves = m_bSplineCurves.size();

    for (int i = 0; i < numberOfSurfaces; i++)
        delete m_bSplineSurfaces[i];

    for (int i = 0; i < numberOfCurves; i++)
        delete m_bSplineCurves[i];

}

void BSplinePatchNetwork::getBoundingBox(int &minX, int &maxX, int &minY, int &maxY, int &minZ, int &maxZ)
{
    const int numCtrPoints = m_properties->getNumberOfControlPoints();
    const QVector3D dummyPoint = m_bSplineSurfaces[0]->getControlPointsP()->at(0).at(0);
    minX = dummyPoint.x();
    maxX = minX;
    minY = dummyPoint.y();
    maxY = minY;
    minZ = dummyPoint.z();
    maxZ = minZ;
    foreach (BSplineSurface *surf, m_bSplineSurfaces) {
        for (int i = 0; i < numCtrPoints; i++) {
            for (int j = 0; j < numCtrPoints; j++) {
                const QVector3D ctrPoint = surf->getControlPointsP()->at(i).at(j);
                if (ctrPoint.x() < minX)
                    minX = ctrPoint.x();
                if (ctrPoint.x() > maxX)
                    maxX = ctrPoint.x();
                if (ctrPoint.y() < minY)
                    minY = ctrPoint.y();
                if (ctrPoint.y() > maxY)
                    maxY = ctrPoint.y();
                if (ctrPoint.z() < minZ)
                    minZ = ctrPoint.z();
                if (ctrPoint.z() > maxZ)
                    maxZ = ctrPoint.z();
            }
        }
    }
}

void BSplinePatchNetwork::translate(const QVector3D &translationVector)
{
    foreach (BSplineSurface *surf, m_bSplineSurfaces)
        surf->translate(translationVector);

    foreach (BSplineCurve *curv, m_bSplineCurves)
        curv->translate(translationVector);
}

void BSplinePatchNetwork::scale(double factorX, double factorY, double factorZ)
{
    foreach (BSplineSurface *surf, m_bSplineSurfaces)
        surf->scale(factorX, factorY, factorZ);

    foreach (BSplineCurve *curv, m_bSplineCurves)
        curv->scale(factorX, factorY, factorZ);
}

CellMesh *BSplinePatchNetwork::getCellMesh()
{
    return m_cellMesh;
}

int BSplinePatchNetwork::getNumberOfBSplineSurfaces()
{
    return m_bSplineSurfaces.size();
}

int BSplinePatchNetwork::getNumberOfBSplineCurves()
{
    return m_bSplineCurves.size();
}

BSplineProperties *BSplinePatchNetwork::getProperties()
{
    return m_properties;
}

BSplineSurface *BSplinePatchNetwork::getBSplineSurface(int id)
{
    return m_bSplineSurfaces[id];
}

QVector<BSplineSurface *> BSplinePatchNetwork::getBSplineSurfaces()
{
    return m_bSplineSurfaces;
}

BSplineCurve *BSplinePatchNetwork::getBSplineCurve(int id)
{
    return m_bSplineCurves[id];
}

QVector<BSplineCurve *> BSplinePatchNetwork::getBSplineCurves()
{
    return m_bSplineCurves;
}

void BSplinePatchNetwork::setBSplines(QVector<BSplineSurface *> &bSplineSurfaces, QVector<BSplineCurve *> &bSplineCurves)
{
    m_bSplineSurfaces = bSplineSurfaces;
    m_bSplineCurves = bSplineCurves;
}

void BSplinePatchNetwork::setBSplines(QVector<BSplineSurface *> &bSplineSurfaces)
{
    m_bSplineSurfaces = bSplineSurfaces;
    this->extractBorderCurves();
}

void BSplinePatchNetwork::addInfo(QString info)
{
    m_info.append(info);
}

void BSplinePatchNetwork::saveToFile(QString filename, QString optionalHeader)
{
    const int numberOfCurves = m_bSplineCurves.size();
    const int numberOfSurfaces = m_bSplineSurfaces.size();
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    out << m_info;
    if (!m_info.endsWith("\n"))
        out << "\n";
    out << optionalHeader << "\n";
    out << "===\n";
    out << QString::number(numberOfSurfaces) << " "
        << QString::number(numberOfCurves) << " "
        << QString::number(m_properties->getOrder()) << " "
        << QString::number(numberOfControlPointsPerRow) << "\n";
    out << "===\n";

    for (int i = 0; i < numberOfSurfaces; i++) {
        QVector<QVector<QVector3D> > *ctrPoints = m_bSplineSurfaces[i]->getControlPointsP();

        for (int u = 0; u < numberOfControlPointsPerRow; u++) {
            for (int v = 0; v < numberOfControlPointsPerRow; v++) {
                QVector3D ctrPoint = ctrPoints->at(u).at(v);
                out << QString::number(ctrPoint.x()) << " " << QString::number(ctrPoint.y()) << " " << QString::number(ctrPoint.z()) << "\n";
            }
        }
        out << "---\n";
    }
    out << "===\n";

    for (int i = 0; i < numberOfCurves; i++) {
        const QVector<QVector3D> *ctrPoints = m_bSplineCurves[i]->getControlPointsP();

        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
            QVector3D ctrPoint = ctrPoints->at(j);
            out << QString::number(ctrPoint.x()) << " " << QString::number(ctrPoint.y()) << " " << QString::number(ctrPoint.z()) << "\n";
        }
        out << "---\n";
    }
    out << "===\n";

    file.close();
}

BSplinePatchNetwork *BSplinePatchNetwork::loadFromFile(QString filename, CellMesh *cm)
{
    QFile file(filename);
    if (!cm || !file.open(QIODevice::ReadOnly))
        return 0;
    QTextStream in(&file);

    QString info;
    QString line = in.readLine();
    while (line != "===") {
        info.append(line).append("\n");
        line = in.readLine();
    }

    //Read general information
    line = in.readLine();
    QStringList split = line.split(" ");
    const int numberOfSurfaces = split[0].toInt();
    const int numberOfCurves = split[1].toInt();
    const int order = split[2].toInt();
    const int numberOfControlPointsPerRow = split[3].toInt();

    line = in.readLine();   //next separator

    //load surfaces
    QVector<BSplineSurface *> surfaces(numberOfSurfaces);
    BSplineProperties *properties = BSplineProperties::createInstanceUniformKnotVectorMultipleEnds(order,
                                                                                                   numberOfControlPointsPerRow,
                                                                                                   cm->getParameterInterval().first,
                                                                                                   cm->getParameterInterval().second);
    for (int i = 0; i < numberOfSurfaces; i++) {
        QVector<QVector<QVector3D> > controlPoints(numberOfControlPointsPerRow, QVector<QVector3D>(numberOfControlPointsPerRow));
        for (int u = 0; u < numberOfControlPointsPerRow; u++) {
            for (int v = 0; v < numberOfControlPointsPerRow; v++) {
                line = in.readLine();
                QStringList tmp = line.split(" ");
                controlPoints[u][v] = QVector3D(tmp[0].toDouble(), tmp[1].toDouble(), tmp[2].toDouble());
            }
        }
        surfaces[i] = new BSplineSurface(properties, properties, controlPoints);
        line = in.readLine();   //next separator ---
    }

    line = in.readLine();   //separator ===

    //load curves
    QVector<BSplineCurve *> curves(numberOfCurves);
    for (int i = 0; i < numberOfCurves; i++) {
        QVector<QVector3D> controlPoints(numberOfControlPointsPerRow);
        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
             line = in.readLine();
             QStringList tmp = line.split(" ");
             controlPoints[j] = QVector3D(tmp[0].toDouble(), tmp[1].toDouble(), tmp[2].toDouble());
        }

        curves[i] = new BSplineCurve(properties, controlPoints);
        line = in.readLine();   //next separator
    }
    file.close();

    BSplinePatchNetwork *result = new BSplinePatchNetwork(cm, properties, surfaces, curves);
    result->addInfo(info);

    return result;
}

void BSplinePatchNetwork::extractBorderCurves()
{
    const int numberOfExistingCurves = m_bSplineCurves.size();
    for (int i = 0; i < numberOfExistingCurves; i++)
        delete m_bSplineCurves[i];
    m_bSplineCurves.clear();

    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();

    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = m_cellMesh->getEdge(i);
        QVector<QVector3D> controlPoints;

        const int fId = ce->getIncidentFaceId(0);
        const int internalCeId = m_cellMesh->getFace(fId)->getInternalEdgeId(ce);
        const bool edgeIsInverted = m_cellMesh->getFace(fId)->edgeIsInverted(internalCeId);

        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
            int id = j;
            if (edgeIsInverted)
                id = numberOfControlPointsPerRow - 1 - j;
            if (internalCeId == 0) {
                controlPoints << m_bSplineSurfaces.at(fId)->getControlPointsP()->at(id).at(0);
            } else if (internalCeId == 1) {
                controlPoints << m_bSplineSurfaces.at(fId)->getControlPointsP()->at(numberOfControlPointsPerRow - 1).at(id);
            } else if (internalCeId == 2) {
                controlPoints << m_bSplineSurfaces.at(fId)->getControlPointsP()->at(id).at(numberOfControlPointsPerRow - 1);
            } else if (internalCeId == 3) {
                controlPoints << m_bSplineSurfaces.at(fId)->getControlPointsP()->at(0).at(id);
            }
        }

        m_bSplineCurves.append(new BSplineCurve(m_properties, controlPoints));
    }
}


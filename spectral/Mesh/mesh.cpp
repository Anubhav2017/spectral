#include "mesh.h"

#include "external/Miniball.hpp"

#include <QFile>
#include <QTextStream>
#include <QStringList>

#include <iostream>
#include <fstream>
#include <cmath>

uint qHash(const QVector3D &v)
{
    return qHash( QString( "%1x%2x%3" ).arg(v.x()).arg(v.y()).arg(v.z()) ) ;
}

Mesh::Mesh()
{

}

Mesh::~Mesh()
{
    foreach (Vertex *v, m_vertices)
        delete v;
    foreach (Edge *e, m_edges)
        delete e;
    foreach (Face *f, m_faces)
        delete f;
}

Face *Mesh::addFace(QVector<QVector3D> corners, QVector3D normal)
{
    const int numberOfCorners = corners.size();

    //Check if vertices exist
    QVector<Vertex *> vertices(numberOfCorners, 0);
    for (int i = 0; i < numberOfCorners; i++) {
        if (m_coordinateVertexMap.contains(corners[i])) {
            vertices[i] = m_coordinateVertexMap[corners[i]];
        } else {
            vertices[i] = new Vertex(corners[i], m_vertices.size());
            m_vertices.append(vertices[i]);
            m_coordinateVertexMap[corners[i]] = vertices[i];
        }
    }

    //Check if corners exist
    QVector<Edge *> edges (numberOfCorners,0);
    for (int i = 0; i < numberOfCorners; i++) {
        int next = (i + 1) % numberOfCorners;
        if (m_vertexPairEdgeMap.contains(QPair<Vertex *, Vertex *>(vertices[i], vertices[next]))) {
            edges[i] = m_vertexPairEdgeMap[QPair<Vertex *, Vertex *>(vertices[i], vertices[next])];
        } else if (m_vertexPairEdgeMap.contains(QPair<Vertex *, Vertex *>(vertices[next], vertices[i]))) {
            edges[i] = m_vertexPairEdgeMap[QPair<Vertex *, Vertex *>(vertices[next], vertices[i])];
        } else {
            int newEdgeId = m_edges.size();
            edges[i] = new Edge(vertices[i], vertices[next], newEdgeId);
            vertices[i]->addIncidentEdgeId(newEdgeId);
            vertices[next]->addIncidentEdgeId(newEdgeId);
            m_edges.append(edges[i]);
            m_vertexPairEdgeMap[QPair<Vertex *, Vertex *>(vertices[i], vertices[next])] = edges[i];
        }
    }

    int newFaceId = m_faces.size();
    Face *fNew = new Face(edges, normal, newFaceId);
    m_faces.append(fNew);
    for (int i = 0; i < numberOfCorners; i++) {
        edges[i]->addIncidentFaceId(newFaceId);
        vertices[i]->addIncidentFaceId(newFaceId);
    }

    return fNew;
}

int Mesh::getNumberOfVertices()
{
    return m_vertices.size();
}

Vertex *Mesh::getVertex(int id)
{
    return m_vertices[id];
}

int Mesh::getNumberOfEdges()
{
    return m_edges.size();
}

Edge *Mesh::getEdge(int id)
{
    return m_edges[id];
}

int Mesh::getNumberOfFaces()
{
    return m_faces.size();
}

Face *Mesh::getFace(int id)
{
    return m_faces[id];
}

Face *Mesh::getNeighbor(Face *f, int id)
{
    Edge * e = f->getEdge(id);
    int id0 = e->getIncidentFaceId(0);
    int id1 = e->getIncidentFaceId(1);
    if (id0 == -1 || id1 == -1)
        return 0;
    else if (id0 == f->getId())
        return m_faces[id1];
    else if (id1 == f->getId())
        return m_faces[id0];
    else
        return 0;  //If the mesh datastructure is used correctly, this should never happen
}

QVector3D Mesh::getVertexNormal(int vertexId)
{
    Vertex *v = m_vertices[vertexId];
    const int numberOfIncidentFaces = v->getIncidentFaceIds()->size();
    QVector3D normal(0, 0, 0);
    for (int i = 0; i < numberOfIncidentFaces; i++)
        normal += m_faces[v->getIncidentFaceIds()->at(i)]->getNormal();

    normal /= (double) numberOfIncidentFaces;
    return normal;
}

QString Mesh::getInfo()
{
    return m_info;
}

void Mesh::resetInfo()
{
    m_info = QString();
}

void Mesh::addInfo(QString additionalInfo)
{
    m_info.append(additionalInfo);
}

int Mesh::getMaxNumberOfEdgesOnVertex()
{
    const int numVertices = m_vertices.size();
    int maxNum = -1;

    for (int i = 0; i < numVertices; i++)
        if (m_vertices[i]->getIncidentEdgeIds()->size() > maxNum)
            maxNum = m_vertices[i]->getIncidentEdgeIds()->size();

    return maxNum;
}

void Mesh::calculateMinMaxVertexPositions(double &minX, double &maxX, double &minY, double &maxY, double &minZ, double &maxZ)
{
    QVector3D pos = m_vertices[0]->getPosition();
    minX = pos.x();
    maxX = minX;
    minY = pos.y();
    maxY = minY;
    minZ = pos.z();
    maxZ = minZ;
    const int numberOfVertices = m_vertices.size();
    for (int i = 1; i < numberOfVertices; i++) {
        pos = m_vertices[i]->getPosition();
        if (pos.x() < minX)
            minX = pos.x();
        if (pos.x() > maxX)
            maxX = pos.x();
        if (pos.y() < minY)
            minY = pos.y();
        if (pos.y() > maxY)
            maxY = pos.y();
        if (pos.z() < minZ)
            minZ = pos.z();
        if (pos.z() > maxZ)
            maxZ = pos.z();
    }
}

void Mesh::calculateMinimumBoundingSphere(double &radius, QVector3D &center)
{
    const int numberOfPoints = m_vertices.size();

    //prepare data
    double **arrayOfPoints = new double*[numberOfPoints];
    for (int i = 0; i < numberOfPoints; i++) {
        double *point = new double[3];
        const QVector3D coordinates = m_vertices[i]->getPosition();
        point[0] = coordinates.x();
        point[1] = coordinates.y();
        point[2] = coordinates.z();
        arrayOfPoints[i] = point;
    }

    typedef double* const* PointIterator;
    typedef const double* CoordIterator;

    //do miniball algorithm
    Miniball::Miniball<Miniball::CoordAccessor<PointIterator, CoordIterator> > minBall(3, arrayOfPoints, arrayOfPoints + numberOfPoints);

    //extract results
    center.setX(minBall.center()[0]);
    center.setY(minBall.center()[1]);
    center.setZ(minBall.center()[2]);

    //Compensate for inaccuracies
    double radiusSquared = minBall.squared_radius();
    for (int i = 0; i < numberOfPoints; i++) {
        const double distSquared = (center - m_vertices[i]->getPosition()).lengthSquared();
        if (distSquared > radiusSquared)
            radiusSquared = distSquared;
    }

    radius = sqrt(radiusSquared);

    //clean up
    for (int i = 0; i < numberOfPoints; i++)
        delete[] arrayOfPoints[i];
    delete[] arrayOfPoints;
}

void Mesh::scaleMeshToFitIntoCube(bool keepProportions, double cubeEdgeLength)
{
    const int numberOfVertices = m_vertices.size();
    double minX, maxX, minY, maxY, minZ, maxZ;
    this->calculateMinMaxVertexPositions(minX, maxX, minY, maxY, minZ, maxZ);

    double difX = maxX - minX;
    double difY = maxY - minY;
    double difZ = maxZ - minZ;
    if (keepProportions) {
        double max = difX;
        if (difY > max)
            max = difY;
        if (difZ > max)
            max = difZ;
        difX = max;
        difY = max;
        difZ = max;
    }

    for (int i = 0; i < numberOfVertices; i++) {
        QVector3D pos = m_vertices[i]->getPosition();

        pos.setX(cubeEdgeLength*(pos.x() - minX)/difX);
        pos.setY(cubeEdgeLength*(pos.y() - minY)/difY);
        pos.setZ(cubeEdgeLength*(pos.z() - minZ)/difZ);

        m_vertices[i]->updatePosition(pos);
    }

    this->addInfo(QString("Scaled to cube with length %1. Keeping ratio: %2\n").arg(cubeEdgeLength).arg(keepProportions));
}

void Mesh::scaleMeshToFitIntoSphere(double radius, QVector3D center)
{
    const int numberOfVertices = m_vertices.size();
    double minBallRadius;
    QVector3D minBallCenter;

    this->calculateMinimumBoundingSphere(minBallRadius, minBallCenter);

    //translate by -minBallCenter, scale my radius/minBallRadius, tanslate by center
    for (int i = 0; i < numberOfVertices; i++) {
        QVector3D pos = m_vertices[i]->getPosition();

        pos.setX((pos.x() - minBallCenter.x()) * radius/minBallRadius + center.x());
        pos.setY((pos.y() - minBallCenter.y()) * radius/minBallRadius + center.y());
        pos.setZ((pos.z() - minBallCenter.z()) * radius/minBallRadius + center.z());

        m_vertices[i]->updatePosition(pos);
    }
}

double Mesh::calculateTotalSurfaceArea()
{
    double result = 0;
    foreach (Face *f, m_faces)
        result += f->calculateAreaByCorners();

    return result;
}

void Mesh::smoothSurface(const int maxIterations, const double maxDistFromOrigin, const double abortTotalDistance)
{
    const int numberOfVertices = m_vertices.size();
    QVector<QVector3D> originalPositions(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++)
        originalPositions[i] = m_vertices[i]->getPosition();

    double totalDist = abortTotalDistance * 2;
    for (int i = 0; i < maxIterations && totalDist > abortTotalDistance; i++) {
        totalDist = 0;
        QVector<QVector3D> newPositions(numberOfVertices);
        for (int i = 0; i < numberOfVertices; i++) {
            Vertex *v = m_vertices[i];
            const int numberOfNeighbors = v->getIncidentEdgeIds()->size();
            QVector3D newPos = v->getPosition()*numberOfNeighbors;
            foreach (int eId, *v->getIncidentEdgeIds())
                newPos += m_edges[eId]->getConnectedVertex(v)->getPosition();

            newPos /= (double) (2*numberOfNeighbors);
            QVector3D translationFromOrigin = newPos - originalPositions[i];
            double distToOrigin = translationFromOrigin.length();
            if (distToOrigin > maxDistFromOrigin) {
                translationFromOrigin.normalize();
                translationFromOrigin *= maxDistFromOrigin;
                newPos = originalPositions[i] + translationFromOrigin;
            }
            totalDist += (newPos - v->getPosition()).length();
            newPositions[i] = newPos;
        }
        totalDist /= numberOfVertices;
        for (int i = 0; i < numberOfVertices; i++)
            m_vertices[i]->updatePosition(newPositions[i]);
    }

    const int numberOfFaces = m_faces.size();
    for (int i = 0; i < numberOfFaces; i++)
        m_faces[i]->recomputeNormal();

    m_info.append(QString("surface smoothed. iterations: %1, maxDistFromOrigin: %2, abortDistance: %3\n")
                  .arg(maxIterations).arg(maxDistFromOrigin).arg(abortTotalDistance));
}

void Mesh::saveToFile(QString filename, QString optionalHeader)
{
    const int numVertices = m_vertices.size();
    const int numEdges = m_edges.size();
    const int numFaces = m_faces.size();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    out << m_info;
    if (!optionalHeader.isEmpty()) {
        out << optionalHeader;
        if (!optionalHeader.endsWith(QString("\n")))
            out << "\n";
    }
    out << "===\n";
    out << QString::number(numVertices) << " " << QString::number(numEdges) << " " << QString::number(numFaces) << "\n";
    out << "===\n";
    for (int i = 0; i < numVertices; i++) {
        Vertex *v = m_vertices[i];
        out << QString::number(v->getId()) << " " <<
               QString::number(v->getPosX()) << " " <<
               QString::number(v->getPosY()) << " " <<
               QString::number(v->getPosZ()) << "\n";
    }

    out << "===\n";

    for (int i = 0; i < numEdges; i++) {
        Edge *e = m_edges[i];
        out << QString::number(e->getId()) << " " <<
               QString::number(e->getVertex(0)->getId()) << " " <<
               QString::number(e->getVertex(1)->getId()) << "\n";
    }

    out << "===\n";

    for (int i = 0; i < numFaces; i++) {
        Face *f = m_faces[i];
        const int numFaceVert = f->getNumberOfVertices();

        out << QString::number(f->getId()) << " " <<
               QString::number(f->getNormalX()) << " " <<
               QString::number(f->getNormalY()) << " " <<
               QString::number(f->getNormalZ()) << " " <<
               QString::number(numFaceVert);
        for (int j = 0; j < numFaceVert; j++)
            out << " " << QString::number(f->getVertex(j)->getId());

        for (int j = 0; j < numFaceVert; j++)
            out << " " << QString::number(f->getEdge(j)->getId());

        out << "\n";
    }

    out << "===";

    file.close();
}

Mesh *Mesh::loadFromFile(QString filename)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly))
        return 0;
    QTextStream in(&file);

    Mesh *m = new Mesh();

    QString info;
    QString line = in.readLine();
    //skip until first seperator
    while (line != "===") {
        info.append(line).append("\n");
        line = in.readLine();
    }
    m->addInfo(info);

    //Read mesh dimensions
    line = in.readLine();   //This line contains the size information of the mesh
    QStringList split = line.split(" ");
    const int numVertices = split[0].toInt();
    const int numEdges = split[1].toInt();
    const int numFaces = split[2].toInt();
    m->m_vertices = QVector<Vertex *>(numVertices, 0);
    m->m_edges = QVector<Edge *>(numEdges, 0);
    m->m_faces = QVector<Face *>(numFaces, 0);

    line = in.readLine();   //next separator

    //read vertices
    line = in.readLine();   //read the first vertex
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        QVector3D pos = QVector3D(split[1].toFloat(), split[2].toFloat(), split[3].toFloat());
        Vertex *v = new Vertex(pos, id);
        m->m_vertices[id] = v;
        m->m_coordinateVertexMap[pos] = v;

        line = in.readLine();
    }

    //read edges
    line = in.readLine();   //read the first edge
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        const int v0Id = split[1].toInt();
        const int v1Id = split[2].toInt();
        Vertex *v0 = m->m_vertices[v0Id];
        Vertex *v1 = m->m_vertices[v1Id];
        v0->addIncidentEdgeId(id);
        v1->addIncidentEdgeId(id);

        Edge *e = new Edge(v0, v1, id);
        m->m_edges[id] = e;
        m->m_vertexPairEdgeMap[QPair<Vertex *, Vertex *>(v0, v1)] = e;

        line = in.readLine();
    }

    //read faces
    line = in.readLine();   //read the first face
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        QVector3D normal(split[1].toFloat(), split[2].toFloat(), split[3].toFloat());
        const int numFaceVert = split[4].toInt();
        QVector<Edge *> faceEdges(numFaceVert);
        for (int i = 0; i < numFaceVert; i++) {
            Edge *e = m->m_edges[split[5 + numFaceVert + i].toInt()];
            faceEdges[i] = e;

            m->m_vertices[split[5 + i].toInt()]->addIncidentFaceId(id); //add incident face to vertex
            e->addIncidentFaceId(id);                                   //add incident face to edge
        }

        m->m_faces[id] = new Face(faceEdges, normal, id);
        line = in.readLine();
    }

    file.close();

    return m;
}

Mesh *Mesh::loadFromBinaryStl(QString filename)
{
    std::ifstream file(filename.toStdString().c_str(), std::ifstream::in | std::ifstream::binary);

    char header_info[80] = "";
    unsigned long nTri = 0;

    //read 80 byte header
    if (file)
        file.read(header_info, 80);
    else
        std::cout << "error" << std::endl;

    //read 4-byte ulong
    if (file)
        file.read((char *) &nTri, 4);
    else
        std::cout << "error" << std::endl;

    //now read in all the triangles
    Mesh *m = new Mesh();
    QString stlheader(header_info);
    m->addInfo("imported from stl. Header: " + stlheader + "\n");
    for(unsigned int i = 0; i < nTri; i++){
        if (file) {
            //read one 50-byte triangle
            float values[12];
            char dummy[2];
            file.read((char *) values, 48);
            file.read(dummy, 2);
            QVector3D normal = QVector3D(values[0], values[1], values[2]);
            QVector<QVector3D> corners(3);
            corners[0] = QVector3D(values[3], values[4], values[5]);
            corners[1] = QVector3D(values[6], values[7], values[8]);
            corners[2] = QVector3D(values[9], values[10], values[11]);

            //add a new triangle to the mesh
            m->addFace(corners, normal);
        }
    }

    file.close();

    return m;
}

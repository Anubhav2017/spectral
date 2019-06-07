#include "plywriter.h"

#include <QFile>
#include <QTextStream>

void PlyWriter::writeBSplineTriangulationWithVertexColors(BSplinePatchNetwork *bSplines, BSplineTriangulation *triangulation, const QVector<QColor> &colormap, const QString filename)
{
    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
        return;

    QTextStream out(&file);

    QString comment("B-Spline triangulation with vertex color");

    const int numVertices = triangulation->getNumberOfVertices();
    if (numVertices != colormap.size())
        return;

    const int numFaces = triangulation->getNumberOfFaces();

    out << "ply\n"
        << "format ascii 1.0\n"
        << "comment" << comment << "\n"
        << "element vertex " << QString::number(numVertices) << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "property uchar alpha\n"
        << "element face " << QString::number(numFaces) << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (int i = 0; i < numVertices; i++) {
        const int surfId = triangulation->getSurfaceIdOfVertex(i);
        const QVector2D param = triangulation->getParameterValueForVertex(i);
        const QVector3D pos = bSplines->getBSplineSurface(surfId)->evaluate(param.x(), param.y());
        const QColor c = colormap[i];
        out << QString::number(pos.x()) << " " << QString::number(pos.y()) << " " << QString::number(pos.z()) << " "
            << QString::number(c.red()) << " " << QString::number(c.green()) << " " << QString::number(c.blue()) << " " << QString::number(c.alpha()) << "\n";
    }

    for (int i = 0; i < numFaces; i++) {
        QVector<int> face = triangulation->getFace(i);
        const int numFaceVert = face.size();
        out << QString::number(numFaceVert);
        for (int j = 0; j < numFaceVert; j++)
            out << " " << QString::number(face[j]);
        out << "\n";
    }

    file.close();
}

void PlyWriter::writeMeshWithVertexColors(Mesh *m, const QVector<QColor> &colormap, const QString filename)
{
    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
        return;

    QTextStream out(&file);

    QString comment("Mesh with vertex color");

    const int numVertices = m->getNumberOfVertices();
    if (numVertices != colormap.size())
        return;

    const int numFaces = m->getNumberOfFaces();

    out << "ply\n"
        << "format ascii 1.0\n"
        << "comment" << comment << "\n"
        << "element vertex " << QString::number(numVertices) << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "property uchar alpha\n"
        << "element face " << QString::number(numFaces) << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (int i = 0; i < numVertices; i++) {
        Vertex *v = m->getVertex(i);
        const QColor c = colormap[i];
        out << QString::number(v->getPosX()) << " " << QString::number(v->getPosY()) << " " << QString::number(v->getPosZ()) << " "
            << QString::number(c.red()) << " " << QString::number(c.green()) << " " << QString::number(c.blue()) << " " << QString::number(c.alpha()) << "\n";
    }

    for (int i = 0; i < numFaces; i++) {
        Face *f = m->getFace(i);
        const int numFaceVert = f->getNumberOfVertices();
        out << QString::number(numFaceVert);
        for (int j = 0; j < numFaceVert; j++)
            out << " " << QString::number(f->getVertex(j)->getId());
        out << "\n";
    }

    file.close();
}

#include<laplacesolver.h>
#include<fstream>

using namespace std;

QVector<double> computeSpectralEigenvector(Mesh *m, int numberOfEigenVector){
    LaplaceSolver ls=LaplaceSolver(m);

    qDebug()<<"ls initiated";
    ls.decompose2();
//    ls.writeVectors();
    qDebug() << "Mesh decomposition done";

    int ef=numberOfEigenVector;

    Eigen::VectorXcd v=ls.getEigenFunction(ef);

    int N=m->getNumberOfVertices();

    QVector<double > vdash(N,0);

    for(int i=0;i<N;i++){
        vdash[i]=v(i).real();
    }

   return vdash;
}

void saveEigenvectorToTxt(Mesh *m, QVector<double> eigenvector, QString filename){


    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
        return;

    QTextStream out(&file);

    const int numVertices = m->getNumberOfVertices();
    const int numFaces = m->getNumberOfFaces();


    for (int i = 0; i < numVertices; i++) {
        Vertex *v = m->getVertex(i);
        const double c = eigenvector[i];
        out << QString::number(v->getPosX()) << " " << QString::number(v->getPosY()) << " " << QString::number(v->getPosZ()) << " "
            << QString::number(c) << "\n";
    }

    file.close();


}

QVector<int> loadMorseSmaleDecomposition(QString filename){

//    QVector<int> faceCellMap(numberOfMeshFaces);

    ifstream readmap;
    readmap.open(filename.toStdString().c_str());
    QVector<int> vertexmap;
    int a;

    while(readmap >> a){
        vertexmap.append(a);
    }
    qDebug() << vertexmap.size();

    return vertexmap;
}


CellMesh *createCellMeshFromMorseSmale(Mesh *m, QVector<int> vertexmap){

    int numberOfMeshFaces=m->getNumberOfFaces();

    QVector<int> faceCellMap(numberOfMeshFaces);



    int N=m->getNumberOfVertices();

    int min,max=2;

    for(int i=0;i<N;i++){

        if(vertexmap[i] < min)min=vertexmap[i];
        if(vertexmap[i] > max)max=vertexmap[i];
    }
    qDebug()<< "min=" << min;
    qDebug() << "max=" << max;
    qDebug() << "Quadmesh code starts here";
    //qDebug() << vertexmap;
    for(int i=0; i<numberOfMeshFaces ;i++){
            int x0 = vertexmap[m->getFace(i)->getVertex(0)->getId()];
            int x1 = vertexmap[m->getFace(i)->getVertex(1)->getId()];
            int x2 = vertexmap[m->getFace(i)->getVertex(2)->getId()];

            if(x0==x1 || x0==x2){
                faceCellMap[i]=x0;
            }
            else if(x1==x0 || x1==x2){
                faceCellMap[i]=x1;
            }
            else faceCellMap[i]=x0;

        }


    for(int i=0; i < numberOfMeshFaces ;i++){
        Edge* e1=m->getFace(i)->getEdge(0);
        Edge* e2=m->getFace(i)->getEdge(1);
        Edge* e3=m->getFace(i)->getEdge(2);
        int f1id=e1->getConnectedFaceId(i);
        int f2id=e2->getConnectedFaceId(i);
        int f3id=e3->getConnectedFaceId(i);
        int x0=faceCellMap[f1id];
        int x1=faceCellMap[f2id];
        int x2=faceCellMap[f3id];

        if(x0==x1 || x0==x2){
            faceCellMap[i]=x0;
        }
        else if(x1==x0 || x1==x2){
            faceCellMap[i]=x1;
        }
        else faceCellMap[i]=x0;

    }

    //qDebug()<< faceCellMap;

//    QVector<QColor> colorMap=ls.generateColorMap2(faceCellMap);
//    qDebug()<< "Color map generated";
//    ls.writeMeshWithVertexColors(colorMap,QString("vector%1s.ply").arg(ef));



//    double r;
//    QVector3D center;
//    m->calculateMinimumBoundingSphere(r, center);

//    for (int i = 0; i < numberOfMeshFaces; i++) {
//        Face *f = m->getFace(i);
//        if (f->getCenter().x() < center.x() && f->getCenter().y() < center.y())
//            faceCellMap[i] = 0;
//        else if (f->getCenter().x() < center.x() && f->getCenter().y() >= center.y())
//            faceCellMap[i] = 1;
//        else if (f->getCenter().x() >= center.x() && f->getCenter().y() < center.y())
//            faceCellMap[i] = 2;
//        else
//            faceCellMap[i] = 3;
//    }

    int numberOfCells = max-min+1;  //this needs to be set to the actual number of cells
    //int numberOfCells = 4;  //this needs to be set to the actual number of cells


    /** Quadmesh code ends here **/

    //Generate cellmesh
    CellMesh *cm = new CellMesh(m, faceCellMap, numberOfCells);
    return cm;
}
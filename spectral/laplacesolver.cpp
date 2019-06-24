#include "laplacesolver.h"
using namespace std;

double LaplaceSolver::getWeights(Edge *e, Vertex *optionalOuterTriangleVertex=0)
{
        double result = 0;
        for (int i = 0; i < 2; i++) {
            const int fId = e->getIncidentFaceId(i);
            if (fId != -1 && !optionalOuterTriangleVertex) {
                Face *f = m_originalMesh->getFace(fId);
                if (f->getNumberOfEdges() != 3)
                    return 1;

                Vertex *v0 = e->getVertex(0);
                Vertex *v1 = e->getVertex(1);
                Vertex *vOut = 0;
                for (int k = 0; !vOut && k < 3; k++) {
                    vOut = f->getVertex(k);
                    if (vOut == v0 || vOut == v1)
                        vOut = 0;
                }
                QVector3D vOut_v0 = v0->getPosition() - vOut->getPosition();
                QVector3D vOut_v1 = v1->getPosition() - vOut->getPosition();
                const double productLength = vOut_v0.length() * vOut_v1.length();
                const double cos = QVector3D::dotProduct(vOut_v0, vOut_v1)/productLength;
                const double sin = QVector3D::crossProduct(vOut_v0, vOut_v1).length()/productLength;
                result += 0.5 * cos/sin; // = cotangent
            } else if (optionalOuterTriangleVertex) {
                Vertex *v0 = e->getVertex(0);
                Vertex *v1 = e->getVertex(1);
                Vertex *vOut = optionalOuterTriangleVertex;

                QVector3D vOut_v0 = v0->getPosition() - vOut->getPosition();
                QVector3D vOut_v1 = v1->getPosition() - vOut->getPosition();
                const double productLength = vOut_v0.length() * vOut_v1.length();
                const double cos = QVector3D::dotProduct(vOut_v0, vOut_v1)/productLength;
                const double sin = QVector3D::crossProduct(vOut_v0, vOut_v1).length()/productLength;
                result += 0.5 * cos/sin; // = cotangent
            }
        }
        return result;
}
void LaplaceSolver::decompose2(){

    //Eigen::SparseMatrix<double> L;
    N=m_originalMesh->getNumberOfVertices();
    Eigen::MatrixXd V(N,3);
    for(int i=0;i<N;i++){
        V(i,0)=double(m_originalMesh->getVertex(i)->getPosX());
        V(i,1)=double(m_originalMesh->getVertex(i)->getPosY());
        V(i,2)=double(m_originalMesh->getVertex(i)->getPosZ());
    }

    int f=m_originalMesh->getNumberOfFaces();

    Eigen::MatrixXi F(f,3);
    for(int i=0;i<f;i++){

        Face* fi=m_originalMesh->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }
    igl::cotmatrix(V,F,Laplacian);

    int nev=80;
    int ncv=200;

    SparseGenMatProd<double> op(Laplacian);
    //GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op,nev,ncv);
    SymEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op,nev,ncv);

    eigs.init();
    int nconv=eigs.compute();

    eigenvalues = eigs.eigenvalues();

    qDebug()<< eigenvalues.rows()<< " " << eigenvalues.cols();
    qDebug()<<"Hi4";
    qDebug()<<eigs.eigenvectors().rows()<<"x"<<eigs.eigenvectors().cols();

        //eigs.eigenvectors(i).resize(N,1);
    functions=eigs.eigenvectors();



}

LaplaceSolver::LaplaceSolver(Mesh* m):tripletList(0){
    m_originalMesh=m;
    N=m->getNumberOfVertices();

    Laplacian=Eigen::SparseMatrix<double>(N,N);
    weights=QVector<QVector<double> >(N,QVector<double>(N,0));
    numberOfEdges=m->getNumberOfEdges();
    edges=QVector<Edge*>(numberOfEdges);
    for(int i=0;i<10;i++){
        m_colorMapColors.append(QColor(i));
    }
    diagonalentries= QVector<double>(N,0);


    for(int k=0;k<numberOfEdges;k++){
        Vertex* v0= m->getEdge(k)->getVertex(0);
        Vertex* v1= m->getEdge(k)->getVertex(1);
        int i=v0->getId();
        int j=v1->getId();
        weights[i][j]=getWeights(m->getEdge(k));

        tripletList.push_back(T(i,j,-weights[i][j]));
        tripletList.push_back(T(j,i,-weights[i][j]));
        diagonalentries[i]+=weights[i][j];
        diagonalentries[j]+=weights[i][j];
    }

    for(int i=0;i<N;i++){
        tripletList.push_back(T(i,i,diagonalentries[i]));
    }
    Laplacian.setFromTriplets(tripletList.begin(),tripletList.end());




}

void LaplaceSolver::decompose(){

    //qDebug()<<"Hi0";
    //Eigen::EigenSolver<Eigen::SparseMatrix<double> > es;
    //qDebug()<<"Hi1"<<endl;

    //cout<<Laplacian<<endl;
    //es.compute(Laplacian);

    qDebug()<<"Hi3";

    int nev=80;
    int ncv=200;

    SparseGenMatProd<double> op(Laplacian);
    //GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op,nev,ncv);
    SymEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op,nev,ncv);

    eigs.init();
    int nconv=eigs.compute();

    eigenvalues = eigs.eigenvalues();

    qDebug()<< eigenvalues.rows()<< " " << eigenvalues.cols();
    qDebug()<<"Hi4";
    qDebug()<<eigs.eigenvectors().rows()<<"x"<<eigs.eigenvectors().cols();

        //eigs.eigenvectors(i).resize(N,1);
    functions=eigs.eigenvectors();
//        qDebug()<<vec.rows()<<"x"<<vec.cols();


  //      functions.append(eigs.eigenvectors(i));
//    int r=eigs.eigenvectors().rows();
//    int c=eigs.eigenvectors().cols();
//    for(int i=0;i<nev;i++){
//        for(int j=0;j<N;j++){
//            functions[i][j]=vec(j,i);
//        }
//    }
    qDebug()<<"Hi5";
//    for(int j=0;j<N;j++){
//        std::complex<double> a= vec(j,3);
//        qDebug() << a.real() << " "<< a.imag();
//    }
    //cout<<functions[0]<<endl;


    //cout<<functions[0]<<endl;
//    for(int i=0;i<80;i++){
//        functions[i]=es.eigenvectors().col(i);
//    }

//    for(int i=0;i<80;i++){
//        std::cout<<eigenvalues[i]<<'\n';
//    }
//    cout<<(Laplacian*eigs.eigenvectors().col(20))-eigs.eigenvectors().col(20)*eigenvalues[20]<<endl;


}

const Eigen::VectorXcd* LaplaceSolver::getEigenValues() const{
    return &eigenvalues;


}

const Eigen::VectorXcd LaplaceSolver::getEigenFunction(int id) const{

    return functions.col(id);
}

//const Eigen::VectorXcd* LaplaceSolver::getEigenFunction(int id) const{
//    return &(functions[id]);
//}

const  std::complex<double> LaplaceSolver::getEigenValue(int id) const{
    return eigenvalues[id];
}

const QVector<QColor> LaplaceSolver::generateColorMap(int id) const{
//   QVector<std::complex<double> >eigenvector=functions.col(id);
   Eigen::VectorXcd eigenvector=functions.col(id);
   QVector<double> realvector(N,0);
   for(int i=0;i<N;i++){
       realvector[i]=eigenvector(i).real();
   }
   const int numSamplePoints = N;

  // QVector<double> metricValues(numSamplePoints);


//   for (int i = 0; i < numSamplePoints; i++) {
//       const int surfId = triangulation->getSurfaceIdOfVertex(i);
//       const QVector2D param = triangulation->getParameterValueForVertex(i);
//       metricValues[i] = metric->evaluatePoint(bSplines.at(surfId), param.x(), param.y());
//   }

   //min/max is not thread safe -> own for loop
   double min = realvector[0];
   double max = realvector[0];
   for (int i = 1; i < numSamplePoints; i++) {
       const double value = realvector[i];
       if (value < min)
           min = value;
       if (value > max)
           max = value;
   }

   qDebug()<< "min="<<min;

    qDebug()<< "max="<<max;

   const double valueRange = max-min ;

   QVector<QColor> colorMap(numSamplePoints);

   for (int i = 0; i < numSamplePoints; i++) {
       double relativeValue = (realvector[i] - min) / valueRange;

       if (relativeValue <= 0.0f) {
           colorMap[i] = m_colorMapColors[0];
       } else if (relativeValue >= 1.0f) {
           colorMap[i] = m_colorMapColors.last();
       } else {
           int interval = 0;
           while (relativeValue > m_colorMapPercentiles[interval + 1]) //value < 1 = colormapPercentiles.last() (due to above if-check)
               interval++;

           const double lower = m_colorMapPercentiles[interval];
           const double upper = m_colorMapPercentiles[interval + 1];
           const double intervalWidth = upper - lower;
           const double lambda = (relativeValue - lower)/intervalWidth;    //(1 - lambda) * lower + (lambda * upper) = relativeValue

           const QVector3D col0(m_colorMapColors[interval].red(), m_colorMapColors[interval].green(), m_colorMapColors[interval].blue());
           const QVector3D col1(m_colorMapColors[interval + 1].red(), m_colorMapColors[interval + 1].green(), m_colorMapColors[interval + 1].blue());
           const QVector3D col = (1 - lambda) * col0 + lambda * col1;

           colorMap[i] = QColor(floor(col.x() + 0.5), floor(col.y() + 0.5), floor(col.z() + 0.5));
       }
   }

   return colorMap;

}

const QVector<QColor> LaplaceSolver::generateColorMap2(int id) const{
   Eigen::VectorXcd eigenvector=functions.col(id);
   QVector<double> realvector(N,0);
   for(int i=0;i<N;i++){
       realvector[i]=eigenvector(i).real();
   }
   const int numSamplePoints = N;

  // QVector<double> metricValues(numSamplePoints);


//   for (int i = 0; i < numSamplePoints; i++) {
//       const int surfId = triangulation->getSurfaceIdOfVertex(i);
//       const QVector2D param = triangulation->getParameterValueForVertex(i);
//       metricValues[i] = metric->evaluatePoint(bSplines.at(surfId), param.x(), param.y());
//   }

   //min/max is not thread safe -> own for loop
   double min = realvector[0];
   double max = realvector[0];
   for (int i = 1; i < numSamplePoints; i++) {
       const double value = realvector[i];
       if (value < min)
           min = value;
       if (value > max)
           max = value;
   }
//   min=-0.039;
//   max=0.028;

   const double valueRange = max-min ;

   QVector<QColor> colorMap(numSamplePoints);

   for (int i = 0; i < numSamplePoints; i++) {
       double relativeValue = 0;
       if(realvector[i]<min){
           relativeValue=0;
       }
       else if (realvector[i]>max){
           relativeValue=1;
       }
       else relativeValue = (realvector[i] - min) / valueRange;

       colorMap[i] = QColor(int(relativeValue*255), int(relativeValue*255), int(relativeValue*255));
       qDebug()<<relativeValue;
       //colorMap[i] = QColor(0, 100, 123);
       }


   return colorMap;

}
const QVector<QColor> LaplaceSolver::generateColorMap2(QVector<std::complex<double> >eigenvector) const{
   //QVector<std::complex<double> >eigenvector=functions[id];
   QVector<double> realvector(N,0);
   for(int i=0;i<N;i++){
       realvector[i]=eigenvector[i].real();
   }
   const int numSamplePoints = N;

  // QVector<double> metricValues(numSamplePoints);


//   for (int i = 0; i < numSamplePoints; i++) {
//       const int surfId = triangulation->getSurfaceIdOfVertex(i);
//       const QVector2D param = triangulation->getParameterValueForVertex(i);
//       metricValues[i] = metric->evaluatePoint(bSplines.at(surfId), param.x(), param.y());
//   }

   //min/max is not thread safe -> own for loop
   double min = realvector[0];
   double max = realvector[0];
   for (int i = 1; i < numSamplePoints; i++) {
       const double value = realvector[i];
       if (value < min)
           min = value;
       if (value > max)
           max = value;
   }
//   min=-0.039;
//   max=0.028;

   const double valueRange = max-min ;

   QVector<QColor> colorMap(numSamplePoints);

   for (int i = 0; i < numSamplePoints; i++) {
       double relativeValue = 0;
       if(realvector[i]<min){
           relativeValue=0;
       }
       else if (realvector[i]>max){
           relativeValue=1;
       }
       else relativeValue = (realvector[i] - min) / valueRange;

       colorMap[i] = QColor(int(relativeValue*255), int(relativeValue*255), int(relativeValue*255));
       //qDebug()<<relativeValue;
       //colorMap[i] = QColor(0, 100, 123);
       }


   return colorMap;

}
void LaplaceSolver::writeMeshWithVertexColors(const QVector<QColor> &colormap, const QString filename)
{
    Mesh *m=m_originalMesh;
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


const void LaplaceSolver::writeMeshWithVector(const Eigen::VectorXcd vec, const QString filename) const
{
    Mesh *m=m_originalMesh;
    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
        return;

    QTextStream out(&file);

    QString comment("Mesh with vertex color");

    const int numVertices = m->getNumberOfVertices();
    const int numFaces = m->getNumberOfFaces();


    for (int i = 0; i < numVertices; i++) {
        Vertex *v = m->getVertex(i);
        const double c = vec(i).real();
        out << QString::number(v->getPosX()) << " " << QString::number(v->getPosY()) << " " << QString::number(v->getPosZ()) << " "
            << QString::number(c) << "\n";
    }

    file.close();
}



//const Eigen::VectorXcd LaplaceSolver::getEigenFunction(int id)const{
//    return functions.col(id);
//}

void LaplaceSolver::writeVectors(){

    for(int i=0;i<80;i++){
        QString filename="file"+QString::number(i);
        QFile file(filename);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
            return;

        QTextStream out(&file);
        for(int j=0;j<N;j++){
            out<<QString::number(functions(i,j).real())<<"\n";
        }
        file.close();
    }
}


const void LaplaceSolver::writeVector(Eigen::VectorXcd v, QString filename) const{
        QFile file(filename);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
            return;

        QTextStream out(&file);
        for(int j=0;j<N;j++){
            out<<QString::number(v(j).real())<<"\n";
        }
        file.close();

}

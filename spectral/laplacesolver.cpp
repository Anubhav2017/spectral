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
LaplaceSolver::LaplaceSolver(Mesh* m){
    m_originalMesh=m;
    N=m->getNumberOfVertices();
    functions=QVector<Eigen::VectorXcd>(80,Eigen::VectorXcd(N));
    Laplacian=MatrixXd(N,N);
    weights=QVector<QVector<double> >(N,QVector<double>(N,0));
    numberOfEdges=m->getNumberOfEdges();
    edges=QVector<Edge*>(numberOfEdges);

    for(int k=0;k<numberOfEdges;k++){
        Vertex* v0= m->getEdge(k)->getVertex(0);
        Vertex* v1= m->getEdge(k)->getVertex(1);
        int i=v0->getId();
        int j=v1->getId();
        weights[i][j]=getWeights(m->getEdge(k));
        Laplacian(i,j)=-weights[i][j];
        Laplacian(i,i)+=weights[i][j];
        Laplacian(j,j)+=weights[i][j];



    }
}

void LaplaceSolver::decompose(){

    Eigen::EigenSolver<MatrixXd> es;
    es.compute(Laplacian);

    eigenvalues = es.eigenvalues();
    for(int i=0;i<80;i++){
        functions[i]=es.eigenvectors().col(i);
    }

    for(int i=0;i<N;i++){
        std::cout<<eigenvalues[i]<<'\n';
    }

}

VectorXcd* LaplaceSolver::getEigenValues(){
    return &eigenvalues;


}

VectorXcd* LaplaceSolver::getEigenFunction(int id){
    return &(functions[id]);
}

std::complex<double> LaplaceSolver::getEigenValue(int id){
    return eigenvalues[id];
}




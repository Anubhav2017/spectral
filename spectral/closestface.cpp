#include "closestface.h"


QVector<int> closestface::loadFaceMap(Mesh *m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);


    for(int i=0;i<r;i++){
//        for(int j=0;j<c;j++){
//            points(i,j)=pointsvector[i][j];
//        }
        points(i,0)=double(pointsvector[i].x());
        points(i,1)=double(pointsvector[i].y());
        points(i,2)=double(pointsvector[i].z());
    }

    int N=m->getNumberOfVertices();
    int Ndash=points.rows();

    Eigen::MatrixXd V(N,3);
    for(int i=0;i<N;i++){
        V(i,0)=double(m->getVertex(i)->getPosX());
        V(i,1)=double(m->getVertex(i)->getPosY());
        V(i,2)=double(m->getVertex(i)->getPosZ());
    }

    int f=m->getNumberOfFaces();

    Eigen::MatrixXi F(f,3);
    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    //qDebug()<< "Hi1";
    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);
    //std::cout <<"I="<<I;
    //qDebug()<<"Ndash="<<Ndash;
    QVector<int> indices(Ndash,0);
    for(int i=0;i<Ndash;i++){
        indices[i]=I(i);
    }

    return indices;
}

QVector<QVector3D> closestface::loadProjections(Mesh *m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);


    for(int i=0;i<r;i++){
//        for(int j=0;j<c;j++){
//            points(i,j)=pointsvector[i][j];
//        }
        points(i,0)=double(pointsvector[i].x());
        points(i,1)=double(pointsvector[i].y());
        points(i,2)=double(pointsvector[i].z());
    }

    int N=m->getNumberOfVertices();
    int Ndash=points.rows();

    Eigen::MatrixXd V(N,3);
    for(int i=0;i<N;i++){
        V(i,0)=double(m->getVertex(i)->getPosX());
        V(i,1)=double(m->getVertex(i)->getPosY());
        V(i,2)=double(m->getVertex(i)->getPosZ());
    }

    int f=m->getNumberOfFaces();

    Eigen::MatrixXi F(f,3);
    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    //qDebug()<< "Hi1";
    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);
    //std::cout <<"I="<<I;
    //qDebug()<<"Ndash="<<Ndash;
    QVector<QVector3D> projections(Ndash,QVector3D(0,0,0));
    for(int i=0;i<Ndash;i++){
        projections[i]=QVector3D(C(i,0),C(i,1),C(i,2));
    }

    return projections;
}

QVector<QVector<double> > closestface::computeBarycentreCoordinates(Mesh* m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){
//        for(int j=0;j<c;j++){
//            points(i,j)=pointsvector[i][j];
//        }
        points(i,0)=double(pointsvector[i].x());
        points(i,1)=double(pointsvector[i].y());
        points(i,2)=double(pointsvector[i].z());
    }

    int N=m->getNumberOfVertices();
    int Ndash=points.rows();

    Eigen::MatrixXd V(N,3);
    for(int i=0;i<N;i++){
        V(i,0)=double(m->getVertex(i)->getPosX());
        V(i,1)=double(m->getVertex(i)->getPosY());
        V(i,2)=double(m->getVertex(i)->getPosZ());
    }

    int f=m->getNumberOfFaces();

    Eigen::MatrixXi F(f,3);
    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);

    Eigen::MatrixXd L;
    qDebug()<< "Hi1";
    Eigen::MatrixXd Va(Ndash,3),Vb(Ndash,3),Vc(Ndash,3);
    for(int i=0;i< Ndash; i++){
        Face* fi=m->getFace(I(i));
        Vertex* v0=fi->getVertex(0);
        Vertex* v1=fi->getVertex(1);
        Vertex* v2=fi->getVertex(2);

        Va(i,0)=double(v0->getPosX());
        Va(i,1)=double(v0->getPosY());
        Va(i,2)=double(v0->getPosZ());

        Vb(i,0)=double(v1->getPosX());
        Vb(i,1)=double(v1->getPosY());
        Vb(i,2)=double(v1->getPosZ());

        Vc(i,0)=double(v2->getPosX());
        Vc(i,1)=double(v2->getPosY());
        Vc(i,2)=double(v2->getPosZ());
    }



    igl::barycentric_coordinates(C,Va,Vb,Vc,L);

    QVector<QVector<double> > L2(Ndash,QVector<double>(3,0));

    for(int i=0;i<Ndash;i++){
        L2[i][0]=L(i,0);
        L2[i][1]=L(i,1);
        L2[i][2]=L(i,2);
    }

    return L2;

}


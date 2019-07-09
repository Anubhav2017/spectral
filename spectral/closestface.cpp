#include "closestface.h"

QVector<int> closestface::loadFaceMap(Mesh* m1, QVector<QVector3D> pointsvector){

    int faceDim=m1->getFace(0)->getNumberOfVertices();

    if(faceDim==4){
        return loadFaceMapQuadMesh(m1,pointsvector);
    }
    else return loadFaceMapTriMesh(m1,pointsvector);
}


QVector<int> closestface::loadFaceMapTriMesh(Mesh *m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){
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
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);

    QVector<int> indices(Ndash,0);
    for(int i=0;i<Ndash;i++){
        indices[i]=I(i);
    }

    return indices;
}


QVector<int> closestface::loadFaceMapQuadMesh(Mesh* mquad, QVector<QVector3D> pointsvector){

    QVector<int> meshmap(0);

    Mesh* m=new Mesh();

    for (int i = 0; i < mquad->getNumberOfFaces(); i++) {
        Face *f = mquad->getFace(i);

        const QVector3D v0 = f->getVertex(0)->getPosition();
        const QVector3D normal=f->getNormal();

        for (int j = 2; j < f->getNumberOfVertices(); j++) {
            QVector3D v1 = f->getVertex(j - 1)->getPosition();
            QVector3D v2 = f->getVertex(j)->getPosition();
            QVector<QVector3D> corners={v0,v1,v2};

            meshmap.push_back(i);
            m->addFace(corners,normal);

        }
    }

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){

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

    int f=m->getNumberOfFaces();;

    Eigen::MatrixXi F(f,3);

    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
;
    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);

    QVector<int> indices(Ndash,0);

    for(int i=0;i<Ndash;i++){
        indices[i]=meshmap[I(i)];
    }

    return indices;

}

QVector<QVector3D> closestface::loadProjections(Mesh *m, QVector<QVector3D> pointsvector){

    int faceDim=m->getFace(0)->getNumberOfVertices();
    if(faceDim==4){
        return loadProjectionsQuadMesh(m,pointsvector);
    }
    else return loadProjectionsTriMesh(m,pointsvector);

}

QVector<QVector3D> closestface::loadProjectionsTriMesh(Mesh *m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){
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
        for(int j=0;j<3;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);

    QVector<QVector3D> projections(Ndash,QVector3D(0,0,0));

    for(int i=0;i<Ndash;i++){
        projections[i]=QVector3D(C(i,0),C(i,1),C(i,2));
    }

    return projections;
}

QVector<QVector3D> closestface::loadProjectionsQuadMesh(Mesh *mquad, QVector<QVector3D> pointsvector){

    QVector<int> meshmap(0);

    Mesh* m=new Mesh();

    for (int i = 0; i < mquad->getNumberOfFaces(); i++) {
        Face *f = mquad->getFace(i);

        const QVector3D v0 = f->getVertex(0)->getPosition();
        const QVector3D normal=f->getNormal();

        for (int j = 2; j < f->getNumberOfVertices(); j++) {
            QVector3D v1 = f->getVertex(j - 1)->getPosition();
            QVector3D v2 = f->getVertex(j)->getPosition();
            QVector<QVector3D> corners={v0,v1,v2};
//            QVector<Edge *> edges(3);
//            edges[0]= new Edge(v0,v1,0);
//            edges[1]= new Edge(v1,v2,1);
//            edges[2]= new Edge(v2,v0,2);

//            Face* fnew =new Face(edges,normal,track);

            meshmap.push_back(i);
            m->addFace(corners,normal);

            }

        }

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);


    for(int i=0;i<r;i++){

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

    QVector<QVector3D> projections(Ndash,QVector3D(0,0,0));
    for(int i=0;i<Ndash;i++){
        projections[i]=QVector3D(C(i,0),C(i,1),C(i,2));
    }

    return projections;

}


QVector<QVector<double> > closestface::computeBarycentreCoordinates(Mesh* m, QVector<QVector3D> pointsvector){
    int faceDim=m->getFace(0)->getNumberOfVertices();
    if(faceDim==4){
        return computeBarycentreCoordinatesQuadMesh(m,pointsvector);
    }
    else return computeBarycentreCoordinatesTriMesh(m,pointsvector);
}



QVector<QVector<double> > closestface::computeBarycentreCoordinatesTriMesh(Mesh* m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){

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

    int faceDim=m->getFace(0)->getNumberOfVertices();

    Eigen::MatrixXi F(f,faceDim);
    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<faceDim;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;

//     bool success = igl::readOBJ("/u/a/agarwala/Desktop/asTorus.obj",V,F);

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

QVector<QVector<double> > closestface::computeBarycentreCoordinatesQuadMesh(Mesh* m, QVector<QVector3D> pointsvector){

    int r=pointsvector.size();
    int c=3;

    Eigen::MatrixXd points(r,c);

    for(int i=0;i<r;i++){

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

    int faceDim=m->getFace(0)->getNumberOfVertices();

    Eigen::MatrixXi F(f,faceDim);
    for(int i=0;i<f;i++){
        Face* fi=m->getFace(i);
        //qDebug() << fi->getNumberOfVertices();
        for(int j=0;j<faceDim;j++){
            F(i,j)=fi->getVertex(j)->getId();
        }
    }

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(points,V,F,sqrD,I,C);

    Eigen::MatrixXd L;
    qDebug()<< "Hi1";
    Eigen::MatrixXd Va(Ndash,3),Vb(Ndash,3),Vc(Ndash,3),Vd(Ndash,3);
    for(int i=0;i< Ndash; i++){
        Face* fi=m->getFace(I(i));
        Vertex* v0=fi->getVertex(0);
        Vertex* v1=fi->getVertex(1);
        Vertex* v2=fi->getVertex(2);
        Vertex* v3=fi->getVertex(3);

        Va(i,0)=double(v0->getPosX());
        Va(i,1)=double(v0->getPosY());
        Va(i,2)=double(v0->getPosZ());

        Vb(i,0)=double(v1->getPosX());
        Vb(i,1)=double(v1->getPosY());
        Vb(i,2)=double(v1->getPosZ());

        Vc(i,0)=double(v2->getPosX());
        Vc(i,1)=double(v2->getPosY());
        Vc(i,2)=double(v2->getPosZ());

        Vd(i,0)=double(v3->getPosX());
        Vd(i,1)=double(v3->getPosY());
        Vd(i,2)=double(v3->getPosZ());
    }



    igl::barycentric_coordinates(C,Va,Vb,Vc,Vd,L);

    QVector<QVector<double> > L2(Ndash,QVector<double>(4,0));

    for(int i=0;i<Ndash;i++){
        L2[i][0]=L(i,0);
        L2[i][1]=L(i,1);
        L2[i][2]=L(i,2);
        L2[i][3]=L(i,3);
    }

    return L2;

}

QVector<QVector2D> closestface::parametrizationQuadMesh(Mesh* mquad, QVector<QVector3D> pointsvector){


    QVector<int> meshmap(0);

    Mesh* m=new Mesh();

    for (int i = 0; i < mquad->getNumberOfFaces(); i++) {
        Face *f = mquad->getFace(i);

        const QVector3D v0 = f->getVertex(0)->getPosition();
        const QVector3D normal=f->getNormal();

        for (int j = 2; j < f->getNumberOfVertices(); j++) {
            QVector3D v1 = f->getVertex(j - 1)->getPosition();
            QVector3D v2 = f->getVertex(j)->getPosition();
            QVector<QVector3D> corners={v0,v1,v2};
//            QVector<Edge *> edges(3);
//            edges[0]= new Edge(v0,v1,0);
//            edges[1]= new Edge(v1,v2,1);
//            edges[2]= new Edge(v2,v0,2);

//            Face* fnew =new Face(edges,normal,track);

            meshmap.push_back(i);
            m->addFace(corners,normal);

            }

        }

    qDebug()<<"m constructed";

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


    Eigen::MatrixXd L;
    qDebug()<< "Hi1";
    Eigen::MatrixXd Va(Ndash,3),Vb(Ndash,3),Vc(Ndash,3),Vd(Ndash,3);
    for(int i=0;i< Ndash; i++){
        Face* fi=mquad->getFace(meshmap[I(i)]);
        Vertex* v0=fi->getVertex(0);
        Vertex* v1=fi->getVertex(1);
        Vertex* v2=fi->getVertex(2);
        Vertex* v3=fi->getVertex(3);

        Va(i,0)=double(v0->getPosX());
        Va(i,1)=double(v0->getPosY());
        Va(i,2)=double(v0->getPosZ());

        Vb(i,0)=double(v1->getPosX());
        Vb(i,1)=double(v1->getPosY());
        Vb(i,2)=double(v1->getPosZ());

        Vc(i,0)=double(v2->getPosX());
        Vc(i,1)=double(v2->getPosY());
        Vc(i,2)=double(v2->getPosZ());

        Vd(i,0)=double(v3->getPosX());
        Vd(i,1)=double(v3->getPosY());
        Vd(i,2)=double(v3->getPosZ());
    }
    qDebug()<< "Hi2";

    igl::barycentric_coordinates(C,Va,Vb,Vc,Vd,L);

    qDebug()<< "Hi3";

    QVector<QVector2D > parameterVector(Ndash);

    Eigen::MatrixXd parameterGenerator(4,2);
    parameterGenerator<< 0,0,0,1,1,1,1,0;

    Eigen::MatrixXd parameters= L*parameterGenerator;


    for(int i=0;i<Ndash;i++){
        parameterVector[i].setX(parameters(i,0));
        parameterVector[i].setY(parameters(i,1));

    }




    for(int i=0;i<Ndash;i++){
        if(parameterVector[i].x()<0){
            //qDebug()<< parameterVector[i].x();
            parameterVector[i].setX(0);
        }
        if(parameterVector[i].x()  >1){
            //qDebug()<< parameterVector[i].x()-1;
            parameterVector[i].setX(1);
        }

        if(parameterVector[i].y()<0){
            //qDebug()<< parameterVector[i].y();
            parameterVector[i].setY(0);
        }

        if(parameterVector[i].y()>1){
            //qDebug()<< parameterVector[i].y()-1;
            parameterVector[i].setY(1);
        }


    }

    return parameterVector;

}

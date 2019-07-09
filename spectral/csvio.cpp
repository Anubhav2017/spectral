#include "csvio.h"


void csvio::createCsvFromVector(QVector<QVector3D> list, QString filename){
    std::ofstream file;
    file.open (QString(filename+".csv").toStdString());

    for(int i=0;i<list.size();i++){
        file<< list[i].x()<<","<< list[i].y()<<","<< list[i].z()<<"\n";
    }

    file.close();
}

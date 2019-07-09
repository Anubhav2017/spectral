#ifndef CSVIO_H
#define CSVIO_H

#include<QVector>
#include<QVector3D>
#include<fstream>


class csvio
{
public:
    static void createCsvFromVector(QVector<QVector3D> list, QString filename);
};

#endif // CSVIO_H

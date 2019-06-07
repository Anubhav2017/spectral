#ifndef VOXELIZER_H
#define VOXELIZER_H

#include <QString>

#include "BSpline/bsplinepatchnetwork.h"

class Voxelizer
{
public:
    static void splinesToVoxelImage(BSplinePatchNetwork *splineNetwork, int sizeX, int sizeY, int sizeZ, QString filename, int samplePoints = 101);
};

#endif // VOXELIZER_H

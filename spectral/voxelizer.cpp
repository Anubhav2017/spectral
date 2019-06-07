#include "voxelizer.h"
#include "external/IASS.h"

#include <cmath>
#include <QQueue>

struct IntPoint {
    IntPoint(int px, int py, int pz) :
        x(px), y(py), z(pz)
    {}

    int x;
    int y;
    int z;
};

void Voxelizer::splinesToVoxelImage(BSplinePatchNetwork *splineNetwork, int sizeX, int sizeY, int sizeZ, QString filename, int samplePoints)
{
    //initialize data
    const uint8_t dummyValue = 2;
    const uint8_t backgroundValue = 0;
    const uint8_t surfaceValue = 1;
    const uint8_t innerValue = 1;
    uint8_t *data = new uint8_t[sizeX*sizeY*sizeZ];
    for (int x = 0; x < sizeX; x++)
        for (int y = 0; y < sizeY; y++)
            for (int z = 0; z < sizeZ; z++)
                data[z * sizeX * sizeY + y * sizeX + x] = dummyValue;


    //Sample surfaces
    const int numberOfPatches = splineNetwork->getNumberOfBSplineSurfaces();
    for (int iSurf = 0; iSurf < numberOfPatches; iSurf++) {
        BSplineSurface *surf = splineNetwork->getBSplineSurface(iSurf);
        for (int iU = 0; iU < samplePoints; iU++) {
            for (int iV = 0; iV < samplePoints; iV++) {
                const double u = (double) iU / (samplePoints - 1);
                const double v = (double) iV / (samplePoints - 1);

                const QVector3D surfacePoint = surf->evaluate(u, v);
                const int x = floor(surfacePoint.x() + 0.5);
                const int y = floor(surfacePoint.y() + 0.5);
                const int z = floor(surfacePoint.z() + 0.5);

                data[z * sizeX * sizeY + y * sizeX + x] = surfaceValue;
            }
        }
    }

    //Do floodfill for the outside (assuming, that the background is one connected component)
    QQueue<IntPoint> queue;
    queue.enqueue(IntPoint(0, 0, 0));
    data[0] = backgroundValue;

    while (!queue.isEmpty()) {
        IntPoint point = queue.dequeue();
        const int x = point.x;
        const int y = point.y;
        const int z = point.z;

        const int neighborX0 = z * sizeX * sizeY + y * sizeX + (x - 1);
        const int neighborX1 = z * sizeX * sizeY + y * sizeX + (x + 1);
        const int neighborY0 = z * sizeX * sizeY + (y - 1) * sizeX + x;
        const int neighborY1 = z * sizeX * sizeY + (y + 1) * sizeX + x;
        const int neighborZ0 = (z - 1) * sizeX * sizeY + y * sizeX + x;
        const int neighborZ1 = (z + 1) * sizeX * sizeY + y * sizeX + x;

        if (x > 1 && data[neighborX0] == dummyValue) {
            queue.enqueue(IntPoint(x - 1, y, z));
            data[neighborX0] = backgroundValue;
        }

        if (x < sizeX - 1 && data[neighborX1] == dummyValue) {
            queue.enqueue(IntPoint(x + 1, y, z));
            data[neighborX1] = backgroundValue;
        }

        if (y > 1 && data[neighborY0] == dummyValue) {
            queue.enqueue(IntPoint(x, y - 1, z));
            data[neighborY0] = backgroundValue;
        }

        if (y < sizeY - 1 && data[neighborY1] == dummyValue) {
            queue.enqueue(IntPoint(x, y + 1, z));
            data[neighborY1] = backgroundValue;
        }

        if (z > 1 && data[neighborZ0] == dummyValue) {
            queue.enqueue(IntPoint(x, y, z - 1));
            data[neighborZ0] = backgroundValue;
        }

        if (z < sizeZ - 1 && data[neighborZ1] == dummyValue) {
            queue.enqueue(IntPoint(x, y, z + 1));
            data[neighborZ1] = backgroundValue;
        }
    }


    //Remaining dummy values are on the inside of the object
    for (int i = 0; i < sizeX * sizeY * sizeZ; i++)
        if (data[i] == dummyValue)
            data[i] = innerValue;

    //setup header
    const int pixelSize = 1;
    IASS::StHeader header;
    header.type = IASS::MONO;
    header.bpp = 1;
    header.stride.x = 1;
    header.stride.y = sizeX;
    header.stride.z = sizeX * sizeY;
    header.size.x = sizeX;
    header.size.y = sizeY;
    header.size.z = sizeZ;
    header.spacing.x = pixelSize;
    header.spacing.y = pixelSize;
    header.spacing.z = pixelSize;

    QByteArray ba = filename.toLatin1();
    const char *filenameCharP = ba.data();

    //Save data
    if (!IASS::Save(header, data, filenameCharP))
        std::cout << "Error: Writing IASS failed" << std::endl;

    delete data;

    std::cout << "Voxel image writen!" << std::endl;
}

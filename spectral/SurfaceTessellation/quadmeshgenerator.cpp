#include "quadmeshgenerator.h"

#include "blossomgraphmatching.h"

#include <QQueue>
#include <QString>

#include <QDebug>

QuadmeshGenerator::QuadmeshGenerator(Mesh *m):
    m_mesh(m), m_parameterizationParametersSet(false)
{

}

CellMesh *QuadmeshGenerator::computeQuadmeshByDelaunayMatching(bool instantTopologyCheckForVoronoi, int minNumberOfFinalCells, int minNumberOfVoronoiCells)
{
    CellMesh *cm = 0;

    //Compute initial Voronoi tesselation
    SurfaceVoronoiGenerator voronoiGenerator(m_mesh, true);
    while (!cm) {
        while (!cm) {
            voronoiGenerator.computeDiskLikeVoronoiTessellation(minNumberOfVoronoiCells, instantTopologyCheckForVoronoi);
            cm = voronoiGenerator.generateCellMesh(true);
            this->applyParametersToCellMesh(cm);

            Face *newSource = 0;

            //check that at most 3 cells meet in a vertex
            CellVertex *cvUnequal3 = this->findCellVertexWithValencyUnequal3(cm);
            if (cvUnequal3) {
                newSource = this->findFreeVoronoiSourceAroundCellVertex(cm, cvUnequal3, &voronoiGenerator);
                if (!newSource) //technically possible, but should not happen in practice
                    return 0;
            }

            //check that cells only meet on one edge (iterate over all cell edges and find duplicates in incident face ids)
            CellEdge *doubleCe = this->findDoubleCellEdgeBetweenCellFaces(cm);
            if (doubleCe) {
                newSource = this->findFreeVoronoiSourceAlongCellEdge(doubleCe, &voronoiGenerator);
                if (!newSource) //technically possible, but should not happen in practice
                    return 0;
            }

            if (newSource) {
                voronoiGenerator.addSource(newSource);
                delete cm;
                cm = 0;
            }

            //check that dual mesh will contain enough cells
            if (cm) {
                const int numberOfDualCells = cm->getNumberOfVertices();
                if (   numberOfDualCells/2 < minNumberOfFinalCells  //note: that needs to be changed in the future if we also allow vertices of valency != 3
                    || numberOfDualCells % 2 == 1)
                {
                    voronoiGenerator.addSourceFromLargestCellWithHigestDistance();
                    delete cm;
                    cm = 0;
                }
            }
        }

        //create dual mesh
        cm->transformToDualMesh();

        //TODO since in the dual mesh cell vertices and cell faces are swapped, we can technically compute (and check!) the matching before the dualization process
        //Compute matching of triangular cells into quads
        const int numberOfCellFaces = cm->getNumberOfFaces();
        const int numberOfCellEdges = cm->getNumberOfEdges();
        QVector<QPair<int, int> > graphEdges(numberOfCellEdges);
        for(int i = 0; i < numberOfCellEdges; i++) {
            CellEdge *ce = cm->getEdge(i);
            graphEdges[i] = QPair<int, int>(ce->getIncidentFaceId(0), ce->getIncidentFaceId(1));
        }

        QVector<int> matching = BlossomGraphMatching::doBlossomGraphMatching(numberOfCellFaces, graphEdges);

        if ((matching.size() == numberOfCellFaces/2)) {
            //merge triangular cells into quads
            for (int i = 0; i < matching.size(); i++)
                cm->mergeCellsAroundEdge(cm->getEdge(matching[i]));

            cm->cleanUpMeshIndices();

            //merging can lead to degenerate cell vertices again -> remove them and merge incident cells
            const int numberOfCellVertices = cm->getNumberOfVertices();
            int cvId = 0;
            while (cvId < numberOfCellVertices) {
                CellVertex *cv = cm->getVertex(cvId);
                if (cv->getIncidentEdgeIds()->size() == 2) {
                    cm->mergeCellsAroundVertex(cv);
                    cm->cleanUpMeshIndices();
                    cvId = -1;  //will get set to zero before loop starts over
                }
                cvId++;
            }

            //Removal of degenerate vertices can reduce the number of cells again
            //Merging too many degenerate cells may also mess up the cell topology
            if (cm->getNumberOfFaces() < minNumberOfFinalCells ||
                cm->getValencyHistogramCellFaces().at(4) != cm->getNumberOfFaces())
            {
                delete cm;
                cm = 0;
                voronoiGenerator.addSourceFromLargestCellWithHigestDistance();
            }
        } else {
            delete cm;
            cm = 0;
            voronoiGenerator.addSourceFromLargestCellWithHigestDistance();
        }
    }

    cm->addInfo(QString("Computed as pairing of Delaunay base complex. Number of cells: %1, minVoronoiCells: %2, minFinalCells: %3, instantTopologyCheck: %4\n")
                .arg(cm->getNumberOfFaces()).arg(minNumberOfVoronoiCells).arg(minNumberOfFinalCells).arg(instantTopologyCheckForVoronoi));

    cm->straightenAllCellEdges();

    return cm;
}

CellMesh *QuadmeshGenerator::computeQuadmeshBySubdividingVoronoi(bool instantTopologyCheckForVoronoi, bool skipRefinement, int minNumberOfFinalCells, int minNumberOfVoronoiCells)
{
    CellMesh *cm = 0;

    SurfaceVoronoiGenerator voronoiGenerator(m_mesh, true);
    while (!cm) {
        voronoiGenerator.computeDiskLikeVoronoiTessellation(minNumberOfVoronoiCells, instantTopologyCheckForVoronoi);
        cm = voronoiGenerator.generateCellMesh(false);
        this->applyParametersToCellMesh(cm);

        //Check that each cell edge consists of at least two mesh edges (i.e. has at least 3 PLVs)
        const int numberOfCellEdges = cm->getNumberOfEdges();
        for (int i = 0; i < numberOfCellEdges && cm; i++) {
            if (cm->getEdge(i)->getPolylineVertices()->size() == 2) {
                voronoiGenerator.addSource(m_mesh->getFace(cm->getEdge(i)->getPolylineVertices()->at(0)->getSharedFaceId(cm->getEdge(i)->getPolylineVertices()->at(1))));
                delete cm;
                cm = 0;
            }
        }

        if (cm) {
            const int numberOfInitialCells = cm->getNumberOfFaces();
            int expectedFinalNumberOfCells = 0;

            for (int i = 0; i < numberOfInitialCells && expectedFinalNumberOfCells >= 0; i++) {
                const int numEdges = cm->getFace(i)->getNumberOfEdges();

                //Assuming number of subdivisions = 1 so far
                if (numEdges < 2) {   //
                    expectedFinalNumberOfCells = -1;
                    voronoiGenerator.addSourceWithLargestDistanceInCell(i);
                    delete cm;
                    cm = 0;
                } else if (numEdges == 2) {
                       expectedFinalNumberOfCells += 2 * 3;
                } else if (numEdges > 4 && !skipRefinement) {  //if refinement is skiped, the else case will also be true here
                     expectedFinalNumberOfCells += (numEdges - 2)/2 * 4 + (numEdges % 2) * 3;    //combination of quads and triangles
                } else {  //triangles, quads and potentially higher order polygons
                      expectedFinalNumberOfCells += 1 * numEdges;
                }
            }

            if (expectedFinalNumberOfCells < minNumberOfFinalCells) {
                voronoiGenerator.addSourceFromLargestCellWithHigestDistance();
                delete cm;
                cm = 0;
            }
        }
    }

    //Refine cell mesh until only faces with 3 or 4 vertices exist

    //Phase one: Find faces with only one or two neighbors and subdivide them
    bool cellsOk = false;
    while (!cellsOk) {
        cellsOk = true;

        const int numberOfCells = cm->getNumberOfFaces();
        for (int i = 0; i < numberOfCells && cellsOk; i++) {
            CellFace *cf = cm->getFace(i);
            if (cf->getNumberOfEdges() == 1) {
                //Subdivide cell
                cellsOk = false;
                QVector<PolylineVertex *> *edgeVertices = ((CellEdge *) cf->getEdge(0))->getPolylineVertices();

                cm->subdivideCellByCut(i, edgeVertices->at(0), edgeVertices->at(edgeVertices->size()/2));
            } else if (cf->getNumberOfEdges() == 2) {
                //Subdivide cell
                cellsOk = false;
                QVector<PolylineVertex *> *edgeVertices0 = ((CellEdge *) cf->getEdge(0))->getPolylineVertices();
                QVector<PolylineVertex *> *edgeVertices1 = ((CellEdge *) cf->getEdge(1))->getPolylineVertices();

                cm->subdivideCellByCut(i, edgeVertices0->at(edgeVertices0->size()/2), edgeVertices1->at(edgeVertices1->size()/2));
            }
        }
    }

    if (!skipRefinement) {

        //Phase two: Find faces with > 4 vertices and subdivide them
        cellsOk = false;
        while (!cellsOk) {
            cellsOk = true;

            const int numberOfCells = cm->getNumberOfFaces();
            for (int i = 0; i < numberOfCells; i++) {
                CellFace *cf = cm->getFace(i);
                if (cf->getNumberOfEdges() > 4) {
                    //Subdivide cell
                    cellsOk = false;

                    //Find vertex pair with lowest valency
                    const int numVertices = cf->getNumberOfVertices();
                    int lowestIndex = 0;
                    int lowestValency = cf->getVertex(0)->getIncidentEdgeIds()->size() + cf->getVertex(numVertices/2)->getIncidentEdgeIds()->size();
                    for (int j = 1; j < numVertices; j++) {
                        int valency = cf->getVertex(j)->getIncidentEdgeIds()->size() + cf->getVertex((j+numVertices/2) % numVertices)->getIncidentEdgeIds()->size();
                        if (valency < lowestValency) {
                            lowestValency = valency;
                            lowestIndex = j;
                        }
                    }
                    //cut out a quadrilateral cell
                    cm->subdivideCellByCut(i, cf->getPolyVertexCorrespondingToCellVertex((CellVertex *) cf->getVertex(lowestIndex)),
                                              cf->getPolyVertexCorrespondingToCellVertex((CellVertex *) cf->getVertex((lowestIndex + 3)%numVertices)));
                }
            }
        }

        cm->addInfo(QString("Cells with > 4 vertices subdivided. New number of cells: %1\n").arg(cm->getNumberOfFaces()));
    }

    //Phase three: Do central subdivision (neccessary to ensure that mesh is actually a quadmesh)
    cm->doSubdivisionForFullMesh();

    cm->addInfo(QString("Computed by subdividing the Voronoi base complex. Number of cells: %1, minVoronoiCells: %2, minFinalCells: %3, instantTopologyCheck: %4\n")
                .arg(cm->getNumberOfFaces()).arg(minNumberOfVoronoiCells).arg(minNumberOfFinalCells).arg(instantTopologyCheckForVoronoi));
    return cm;
}

CellMesh *QuadmeshGenerator::subdivideMeshByCuttingPlanes(QVector<QVector3D> &normals, QVector<QVector3D> &basePoints)
{
    const int numberOfMeshFaces = m_mesh->getNumberOfFaces();

    const int numberOfCuttingPlanes = normals.size();
    if (numberOfCuttingPlanes != basePoints.size())
        return 0;

    //subdivide mesh along planes
    const int maxNumberOfCells = pow(2, numberOfCuttingPlanes);

    QList<QVector<Face *> > cellsList;
    cellsList.reserve(maxNumberOfCells);
    for(int i = 0; i < maxNumberOfCells; i++)
        cellsList.append(QVector<Face *>());
    QVector<int> faceCellMap(numberOfMeshFaces, -1);

    int cutId[numberOfCuttingPlanes];
    for (int j = 0; j < numberOfCuttingPlanes; j++)
        cutId[j] = pow(2, j);

    for (int i = 0; i < numberOfMeshFaces; i++) {
        int targetId = 0;
        Face *f = m_mesh->getFace(i);
        const QVector3D pos = f->getCenter();
        for (int j = 0; j < numberOfCuttingPlanes; j++)
            if (QVector3D::dotProduct(pos - basePoints[j], normals[j]) > 0)
                targetId += cutId[j];

        cellsList[targetId] << f;
    }

    //clean up empty cells
    for (int j = maxNumberOfCells - 1; j >= 0; j--)
        if (cellsList[j].isEmpty())
            cellsList.removeAt(j);

    //set up face cell map
    const int numberOfCellsAfterComponentSeparation = cellsList.size();
    for (int j = 0; j < numberOfCellsAfterComponentSeparation; j++) {
        const int numFacesInCell = cellsList[j].size();
        for (int i = 0; i < numFacesInCell; i++)
            faceCellMap[cellsList[j][i]->getId()] = j;
    }

    //if a cell is not one connected component, subdivide it into several cells
    const int numberOfCellsAfterCut = cellsList.size();
//    const double infinity = std::numeric_limits<double>::infinity();
    QVector<bool> processed(numberOfMeshFaces, false);
    for (int cellId = 0; cellId < numberOfCellsAfterCut; cellId++) {
        const int numberOfFacesInCell = cellsList.at(cellId).size();
        QList<QVector<Face *> > connectedComponents;

        Face *unprocessedFace = cellsList.at(cellId).first();
        int processedFaces = 0;
        //Do BFS starting at any unprocessed face until all faces have been processed
        while (unprocessedFace) {
            QVector< Face *> component;
            QQueue<Face *> queue;
            queue.enqueue(unprocessedFace);
            processed[unprocessedFace->getId()] = true;
            while (!queue.isEmpty()) {
                Face *f = queue.dequeue();
                component << f;

                const int numberOfNeighbors = f->getNumberOfEdges();
                for (int i = 0; i < numberOfNeighbors; i++) {
                    Face * fNeighbor = m_mesh->getNeighbor(f, i);
                    const int fNeighborId = fNeighbor->getId();
                    if (!processed[fNeighborId] && faceCellMap[fNeighborId] == cellId) {
                        processed[fNeighborId] = true;
                        queue.enqueue(fNeighbor);
                    }
                }
            }

            //find next unprocessed face (if it exists)
            unprocessedFace = 0;
            connectedComponents.append(component);
            processedFaces += component.size();
            if (processedFaces != numberOfFacesInCell)
                for (int i = 0; i < numberOfFacesInCell; i++)
                    if (!processed[cellsList.at(cellId).at(i)->getId()])
                        unprocessedFace = cellsList.at(cellId).at(i);
        }

        //check number of connected components and update global subdivision accordingly
        const int numberOfConnectedComponents = connectedComponents.size();
        if (numberOfConnectedComponents > 1) {
            cellsList[cellId] = connectedComponents[0];

            for (int i = 1; i < numberOfConnectedComponents; i++) {
                int newCellId = cellsList.size();
                cellsList << connectedComponents[i];
                const int numberOfFacesInComponent = connectedComponents[i].size();

                //Update face cell map
                for (int k = 0; k < numberOfFacesInComponent; k++)
                    faceCellMap[connectedComponents[i][k]->getId()] = newCellId;
            }
        }
    }

    //TODO: topology check -> for now just return NULL in case of non-disks
    qDebug() << "TODO Subdivision by cutting plane needs topology check";

    QString basepointsStr = ("Basepoints:");
    foreach (QVector3D basePoint, basePoints)
        basepointsStr.append("\n (" + QString::number(basePoint.x()) + ", " + QString::number(basePoint.y()) + ", " + QString::number(basePoint.z()) + ")");
    QString normalsStr = ("Normals:");
    foreach (QVector3D normal, normals)
        normalsStr.append("\n (" + QString::number(normal.x()) + ", " + QString::number(normal.y()) + ", " + QString::number(normal.z()) + ")");

    CellMesh *cm = new CellMesh(m_mesh, faceCellMap, cellsList.size());
    this->applyParametersToCellMesh(cm);
    cm->addInfo(QString("Computed by cutting planes with\n %1 \n %2 \n resulting in %3 number of cells\n")
                .arg(basepointsStr).arg(normalsStr).arg(cm->getNumberOfFaces()));
    return cm;
}

void QuadmeshGenerator::setParameterizationParameters(Parameterization::BorderType bt, Parameterization::WeightType wt, QPair<double, double> paramInterval, int storageReservation)
{
    m_parameterizationParametersSet = true;
    m_paramBorderType = bt;
    m_paramWeightType = wt;
    m_paramParameterInterval = paramInterval;
    m_paramStorageReservation = storageReservation;
}

CellVertex *QuadmeshGenerator::findCellVertexWithValencyUnequal3(CellMesh *cm)
{
    const int numberOfCellVertices = cm->getNumberOfVertices();
    for (int i = 0; i < numberOfCellVertices; i++) {
        CellVertex *cv = cm->getVertex(i);
        if (cv->getIncidentEdgeIds()->size() != 3)
            return cv;
    }
    return 0;
}

CellEdge *QuadmeshGenerator::findDoubleCellEdgeBetweenCellFaces(CellMesh *cm)
{
    const int numberOfCellEdges = cm->getNumberOfEdges();
    QSet<QPair<int, int> > edgeSet;
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = cm->getEdge(i);
        const int fId0 = ce->getIncidentFaceId(0);
        const int fId1 = ce->getIncidentFaceId(1);
        QPair<int, int> sortedNeighborhoodRelation; //if we have several edges between two faces, fId0 and fId1 might be in different order
        if (fId0 < fId1)
            sortedNeighborhoodRelation = QPair<int, int>(fId0, fId1);
        else if (fId1 < fId0)
            sortedNeighborhoodRelation = QPair<int, int>(fId1, fId0);
        else    //they are equal. Should not happen
            return ce;

        if (!edgeSet.contains(sortedNeighborhoodRelation))
            edgeSet.insert(sortedNeighborhoodRelation);
        else
            return ce;
    }

    return 0;
}

Face *QuadmeshGenerator::findFreeVoronoiSourceAroundCellVertex(CellMesh *cm, CellVertex *cv, SurfaceVoronoiGenerator *voronoiGenerator)
{
    CellEdge *incidentCe = cm->getEdge(cv->getIncidentEdgeIds()->first());
    PolylineVertex *plv = incidentCe->getPolylineVertices()->first();
    if (plv->getPosition() != cv->getPosition())
        plv = incidentCe->getPolylineVertices()->last();

    const int numIncidentFaces = plv->getIncidentFaceIds()->size();
    for (int i = 0; i < numIncidentFaces; i++) {
        Face *f = m_mesh->getFace(plv->getIncidentFaceIds()->at(i));
        if (!voronoiGenerator->isSource(f))
            return f;
    }

    return 0;
}

Face *QuadmeshGenerator::findFreeVoronoiSourceAlongCellEdge(CellEdge *ce, SurfaceVoronoiGenerator *voronoiGenerator)
{
    const int numPLV = ce->getPolylineVertices()->size();
    for (int j = 0; j < numPLV; j++) {
        PolylineVertex *plv = ce->getPolylineVertices()->at((j + numPLV/2) % numPLV);
        const int numIncidentFaces = plv->getIncidentFaceIds()->size();
        for (int i = 0; i < numIncidentFaces; i++) {
            Face *f = m_mesh->getFace(plv->getIncidentFaceIds()->at(i));
            if (!voronoiGenerator->isSource(f))
                return f;
        }
    }

    return 0;
}

void QuadmeshGenerator::applyParametersToCellMesh(CellMesh *cm)
{
    if (m_parameterizationParametersSet) {
        cm->setParameterizationStorageReservation(m_paramStorageReservation);
        cm->setParameterizationWeightType(m_paramWeightType);
        cm->setParameterizationBorderType(m_paramBorderType);
        cm->setParameterInterval(m_paramParameterInterval);
    }
}


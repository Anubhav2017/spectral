#include "cellmesh.h"

#include <cmath>
#include <algorithm>
#include <iostream>

#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QSet>

#include <QDebug>

CellMesh::CellMesh(Mesh *originalMesh, QPair<double, double> parameterInterval):
    m_cellFaces(QVector<CellFace *>(1)), m_cellFaceParameterizations(QVector<Parameterization *>(1, 0)),
    m_parameterInterval(parameterInterval),
    m_originalMesh(originalMesh), m_parameterizationWeightType(Parameterization::WT_cotangent), m_parameterizationBorderType(Parameterization::BT_lengthWeightedCircle),
    m_parameterizationStorageReservation_ExpectedEdgesPerVertex(originalMesh->getMaxNumberOfEdgesOnVertex() + 1),
    m_info(QString("Created with only one cell containing the whole original mesh\n"))
{
    const int numVertices = originalMesh->getNumberOfVertices();
    QVector<Vertex *> meshVertices(numVertices);
    for (int i = 0; i < numVertices; i++)
        meshVertices[i] = originalMesh->getVertex(i);
    const int numEdges = originalMesh->getNumberOfEdges();
    QVector<Edge *> meshEdges(numEdges);
    for (int i = 0; i < numEdges; i++)
        meshEdges[i] = originalMesh->getEdge(i);

    QVector<Edge *> dummyCellEdges;
    CellFace *cf = new CellFace(dummyCellEdges, meshEdges, meshVertices, 0);

    m_cellFaces[0] = cf;
}

CellMesh::CellMesh(Mesh *originalMesh, QVector<int> &faceCellMap, int numberOfCells, bool atLeast3PlvPerEdge, QPair<double, double> parameterInterval):
    m_cellFaces(QVector<CellFace *>(numberOfCells)), m_cellFaceParameterizations(QVector<Parameterization *>(numberOfCells, 0)),
    m_parameterInterval(parameterInterval), m_originalMesh(originalMesh),
    m_parameterizationWeightType(Parameterization::WT_cotangent), m_parameterizationBorderType(Parameterization::BT_lengthWeightedCircle),
    m_parameterizationStorageReservation_ExpectedEdgesPerVertex(originalMesh->getMaxNumberOfEdgesOnVertex() + 1),
    m_info(QString("Created by defining %1 individual cells\n").arg(numberOfCells))
{
    //Iterate once over the mesh vertices and edges and assign them to the correct cells
    QVector<QVector<Vertex *> > containedVertices(numberOfCells, QVector<Vertex *>());
    QVector<QVector<Edge *> > containedEdges(numberOfCells, QVector<Edge *>());

    const int numberOfMeshVertices = originalMesh->getNumberOfVertices();
    QVector<int> vertexCellMap(numberOfMeshVertices);

    //assign mesh vertices to cells
    for (int i = 0; i < numberOfMeshVertices; i++) {
        Vertex *v = originalMesh->getVertex(i);
        const int numIncidentFaces = v->getIncidentFaceIds()->size();
        const int cellId = faceCellMap[v->getIncidentFaceIds()->at(0)];
        bool inside = true;
        for (int j = 1; j < numIncidentFaces && inside; j++)
            inside = (faceCellMap[v->getIncidentFaceIds()->at(j)] == cellId);

        if (inside) {
            vertexCellMap[i] = cellId;
            containedVertices[cellId].append(v);
        } else {
            vertexCellMap[i] = -1;
        }
    }

    //assign mesh edges to cells (or borders)
    QVector<QList<Edge *> > borders(numberOfCells, QList<Edge *>());
    const int numberOfMeshEdges = originalMesh->getNumberOfEdges();
    for (int i = 0; i < numberOfMeshEdges; i++) {
        Edge *e = originalMesh->getEdge(i);
        const int vId = vertexCellMap[e->getVertex(0)->getId()];
        if (vId != -1 && vId == vertexCellMap[e->getVertex(1)->getId()]) {
            containedEdges[vId].append(e);
        } else {
            const int f0Id = faceCellMap[e->getIncidentFaceId(0)];
            const int f1Id = faceCellMap[e->getIncidentFaceId(1)];
            if (f0Id != f1Id) {
                borders[f0Id].append(e);
                borders[f1Id].append(e);
            }
        }
    }

    //sort borders
    for (int i = 0; i < numberOfCells; i++) {
        QList<Edge *> unsortedEdges = borders[i];
        QList<Edge *> sortedEdges;

        Edge *firstEdge = unsortedEdges.takeFirst();
        sortedEdges.append(firstEdge);
        Vertex *startPoint = firstEdge->getVertex(0);
        Vertex *endPoint = firstEdge->getVertex(1);

        while (!unsortedEdges.isEmpty()) {
            int nextIndex = -1;
            Edge *next = 0;
            Vertex *nextEndPoint = 0;
            for (int j = 0; j < unsortedEdges.size(); j++) {
                if (unsortedEdges[j]->getVertex(0) == endPoint) {
                    if (next) {
                        qDebug() << "ERROR in CellMesh construction: border has self intersections";
                        exit(-2);
                    }
                    next = unsortedEdges[j];
                    nextEndPoint = next->getVertex(1);
                    nextIndex = j;
                } else if (unsortedEdges[j]->getVertex(1) == endPoint) {
                    if (next) {
                        qDebug() << "ERROR in CellMesh construction: border has self intersections";
                        exit(-2);
                    }
                    next = unsortedEdges[j];
                    nextEndPoint = next->getVertex(0);
                    nextIndex = j;
                }
            }

            if (nextIndex == -1) {
                //That should not happen -> border is not complete
                qDebug() << "ERROR in CellMesh construction: border of a cell is not connected";
                exit(-2);
            }

            unsortedEdges.removeAt(nextIndex);
            sortedEdges << next;
            endPoint = nextEndPoint;
        }

        if (endPoint != startPoint) {
            qDebug() << "ERROR in CellMesh construction: border of a cell is not closed";
            exit(-2);
        }

        borders[i] = sortedEdges;
    }

    //Now, construct the cell vertices, edges, and faces
    QHash<Vertex *, CellVertex *> vertexCellVertexMap;
    for (int cId = 0; cId < numberOfCells; cId++) {
        //Split border into segments such that each segment separates exactly two cells
        QList<Edge *> *border = &borders[cId];

        QVector<Vertex *> cellVertices;
        QList<QList<Edge *> > borderSegments;

        int neighborId = faceCellMap.at(border->last()->getIncidentFaceId(0));
        if (neighborId == cId)
            neighborId = faceCellMap.at(border->last()->getIncidentFaceId(1));

        QList<Edge *> borderSegment;
        Edge *eLast = border->last();
        int firstVertexId = -1; //Will be set to the occurence of the first vertex in the border -> used to fill the border segment that crosses goes from the last to the first edge
        for (int i = 0; i < border->size(); i++) {
            Edge *e = border->at(i);
            int nextId = faceCellMap.at(e->getIncidentFaceId(0));
            if (nextId == cId)
                nextId = faceCellMap.at(e->getIncidentFaceId(1));

            if (nextId != neighborId) {
                neighborId = nextId;
                Vertex *meshVertex = e->getSharedVertex(eLast);
                CellVertex *v = 0;
                if (vertexCellVertexMap.contains(meshVertex)) {
                    v = vertexCellVertexMap[meshVertex];
                } else {
                    v = new CellVertex(meshVertex->getPosition(), m_cellVertices.size());
                    vertexCellVertexMap[meshVertex] = v;
                    m_cellVertices.append(v);
                }

                cellVertices.append(v);
                borderSegments.append(borderSegment);
                borderSegment = QList<Edge *>();

                if (firstVertexId == -1)
                    firstVertexId = i;
            }

            borderSegment.append(e);
            eLast = e;
        }

        if (borderSegments.isEmpty())
            borderSegments.append(borderSegment);
        else
            for (int i = borderSegment.size()-1; i >= 0; i--)
                borderSegments.first().prepend(borderSegment[i]);

        //Now, create cell edges (if needed)
        const int numberOfNeighbors = borderSegments.size();
        QVector<Edge *> cellEdges(numberOfNeighbors, 0);
        if (numberOfNeighbors > 1) {
            for (int i = 0; i < numberOfNeighbors; i++) {
                Edge *e = borderSegments[i][0];
                int neighborId = faceCellMap.at(e->getIncidentFaceId(0));
                if (neighborId == cId)
                    neighborId = faceCellMap.at(e->getIncidentFaceId(1));

                CellVertex *vStart = (CellVertex *) cellVertices[(cellVertices.size()+i-1)%cellVertices.size()];
                CellVertex *vEnd = (CellVertex *) cellVertices[i];

                CellEdge *cEdge = 0;
                if (neighborId < cId) {
                    //Edge already exists -> take it
                    int j = 0;
                    while (this->getNeighborId(m_cellFaces[neighborId], j) != cId ||
                           !m_cellFaces[neighborId]->getEdge(j)->isIncident(vStart) ||
                           !m_cellFaces[neighborId]->getEdge(j)->isIncident(vEnd)) {
                        j++;
                    }
                    cEdge = (CellEdge *) m_cellFaces[neighborId]->getEdge(j);
                } else {
                    //Create new cell edge
                    int newEdgeId = m_cellEdges.size();
                    QVector<PolylineVertex *> cellEdgeVertices;
                    if (borderSegments[i].size() < 2) {
                        cellEdgeVertices.append(new PolylineVertex(borderSegments[i][0]->getVertex(0)));
                        if (atLeast3PlvPerEdge) {
                            const QVector3D centralPoint = (borderSegments[i][0]->getVertex(0)->getPosition() + borderSegments[i][0]->getVertex(1)->getPosition())/2;
                            cellEdgeVertices.append(new PolylineVertex(centralPoint, borderSegments[i][0]));
                        }
                        cellEdgeVertices.append(new PolylineVertex(borderSegments[i][0]->getVertex(1)));
                    } else {
                        cellEdgeVertices.append(new PolylineVertex(borderSegments[i][0]->getUnsharedVertex(borderSegments[i][1])));
                        const int numEdges = borderSegments[i].size();
                        for (int j = 1; j < numEdges; j++)
                            cellEdgeVertices.append(new PolylineVertex(borderSegments[i][j-1]->getSharedVertex(borderSegments[i][j])));

                        cellEdgeVertices.append(new PolylineVertex(borderSegments[i][numEdges-1]->getUnsharedVertex(borderSegments[i][numEdges-2])));
                    }
                    cEdge = new CellEdge(vStart, vEnd, cellEdgeVertices, newEdgeId);
                    cEdge->addIncidentFaceId(cId);
                    cEdge->addIncidentFaceId(neighborId);
                    vStart->addIncidentEdgeId(newEdgeId);
                    vEnd->addIncidentEdgeId(newEdgeId);
                    m_cellEdges.append(cEdge);
                }

                cellEdges[i] = cEdge;

            }
        } else {
            //Only one neighbor -> no vertices generated by the procedure above
            Edge *e = borderSegments[0][0];
            int neighborId = faceCellMap.at(e->getIncidentFaceId(0));
            if (neighborId == cId)
                neighborId = faceCellMap.at(e->getIncidentFaceId(1));

            CellEdge *cEdge = 0;
            if (neighborId < cId) {
                CellFace *neighbor = m_cellFaces[neighborId];
                int i = 0;
                cEdge = (CellEdge *) neighbor->getEdge(i);
                while (i < neighbor->getNumberOfEdges()-1 && !(cEdge->getIncidentFaceId(0) == cId || cEdge->getIncidentFaceId(1) == cId)) {
                    i++;
                    cEdge = (CellEdge *) neighbor->getEdge(i);
                }
                CellVertex *cv = (CellVertex *) cEdge->getVertex(0);
                cellVertices.append(cv);
            } else {
                Vertex *v = e->getUnsharedVertex(borderSegments[0][1]);

                CellVertex *cv = new CellVertex(v->getPosition(), m_cellVertices.size());
                vertexCellVertexMap[v] = cv;
                m_cellVertices.append(cv);
                cellVertices.append(cv);

                QVector<PolylineVertex *> cellEdgeVertices;
                cellEdgeVertices.append(new PolylineVertex(v));
                const int numEdges = borderSegments[0].size();
                for (int j = 1; j < numEdges; j++)
                    cellEdgeVertices.append(new PolylineVertex(borderSegments[0][j-1]->getUnsharedVertex(borderSegments[0][j])));
                cellEdgeVertices.append(new PolylineVertex(borderSegments[0][numEdges-1]->getUnsharedVertex(borderSegments[0][numEdges-2])));

                int newEdgeId = m_cellEdges.size();
                cEdge = new CellEdge(cv, cv, cellEdgeVertices, newEdgeId);
                cEdge->addIncidentFaceId(cId);
                cEdge->addIncidentFaceId(neighborId);
                cv->addIncidentEdgeId(newEdgeId);
                m_cellEdges.append(cEdge);
            }
            cellEdges[0] = cEdge;
        }

        CellFace *f = new CellFace(cellEdges, containedEdges[cId], containedVertices[cId], cId);
        foreach (Vertex *v, cellVertices)
            v->addIncidentFaceId(cId);

        m_cellFaces[cId] = f;
    }
}

CellMesh::~CellMesh()
{
    const int numberOfCellFaces = m_cellFaces.size();
    const int numberOfCellEdges = m_cellEdges.size();
    const int numberOfCellVertices = m_cellVertices.size();

    for (int i = 0; i < numberOfCellFaces; i++) {
        if (m_cellFaceParameterizations[i])
            delete m_cellFaceParameterizations[i];

        delete m_cellFaces[i];
    }

    for (int i = 0; i < numberOfCellEdges; i++) {
        delete m_cellEdges[i];
    }

    for (int i = 0; i < numberOfCellVertices; i++)
        delete m_cellVertices[i];
}

void CellMesh::straightenAllCellEdges()
{
    this->cleanUpMeshIndices(); //Necessary to avoid NULLs in vertex, edge, and face lists

    const int numberOfEdges = m_cellEdges.size();
    for (int i = 0; i < numberOfEdges; i++)
        this->straightenCellEdgeWithSafetyCheck(m_cellEdges[i]);

    this->addInfo(QString("Edges straightened!\n"));

    this->cleanUpMeshIndices();
}

void CellMesh::optimizeAllVertices()
{
    this->cleanUpMeshIndices(); //Necessary to avoid NULLs in vertex, edge, and face lists

    const int numberOfVertices = m_cellVertices.size();
    for (int i = 0; i < numberOfVertices; i++)
        this->optimizeVertexWithSafetyCheck(m_cellVertices[i]);

    this->addInfo(QString("Cell vertex positions optimized!\n"));

    this->cleanUpMeshIndices();
}

void CellMesh::optimizeAllVerticesWithTemporarySubdivision()
{
    this->cleanUpMeshIndices(); //Necessary to avoid NULLs in vertex, edge, and face lists
    const int numberOfInitialVertices = m_cellVertices.size();
    const int numberOfInitialEdges = m_cellEdges.size();

    //Isolate Cell Vertices
    this->doSubdivisionForFullMesh();

    //Optimize Isolated Vertices
    for (int i = 0; i < numberOfInitialVertices; i++)
        optimizeVertexUnsafe(m_cellVertices[i]);

    //Merge subdivided cells again
    const int newNumberOfCellVertices = m_cellVertices.size();
    for (int i = numberOfInitialVertices + numberOfInitialEdges; i < newNumberOfCellVertices; i++)
        mergeCellsAroundVertex(m_cellVertices[i]);

    //Remove cell vertices along the subdivided edges
    for (int i = numberOfInitialVertices - 1; i < numberOfInitialVertices + numberOfInitialEdges; i++)
        mergeEdgesAroundUnecessaryVertex(m_cellVertices[i]);

    //boundaries are probably not straight anymore, since only the original cell vertices have been optimized
    this->straightenAllCellEdges();

    this->addInfo(QString("Cell vertex positions optimized with temporary additional subdivision!\n"));

    this->cleanUpMeshIndices();
}

CellEdge *CellMesh::straightenCellEdgeWithSafetyCheck(CellEdge *ce)
{
    //check if boundary of surrounding macro cell has self intersections (parameterization cannot handle that)
    bool selfIntersection = false;
    QSet<CellVertex *> cellVerticesOnBoundary;
    CellVertex *cv0 = (CellVertex *) ce->getVertex(0);
    CellVertex *cv1 = (CellVertex *) ce->getVertex(1);
    for (int k = 0; k < 2; k++) {
        CellFace *incidentCf = m_cellFaces[ce->getIncidentFaceId(k)];
        const int numCV = incidentCf->getNumberOfVertices();
        for (int iCv = 0; iCv < numCV; iCv++) {
            CellVertex *incidentCv = (CellVertex *) incidentCf->getVertex(iCv);
            if (incidentCv != cv0 && incidentCv != cv1 && !cellVerticesOnBoundary.contains(incidentCv))
                cellVerticesOnBoundary.insert(incidentCv);
            else if (cellVerticesOnBoundary.contains(incidentCv))
                selfIntersection = true;
        }
    }

    if (selfIntersection)
        return this->straightenCellEdgeSafe(ce);
    else
        return this->straightenCellEdgeUnsafe(ce);
}

CellEdge *CellMesh::straightenCellEdgeUnsafe(CellEdge *ce)
{
    CellVertex *cv0 = (CellVertex *) ce->getVertex(0);
    CellVertex *cv1 = (CellVertex *) ce->getVertex(1);
    CellFace *cf = this->mergeCellsAroundEdge(ce);
    PolylineVertex *plv0 = cf->getPolyVertexCorrespondingToCellVertex(cv0);
    PolylineVertex *plv1 = cf->getPolyVertexCorrespondingToCellVertex(cv1);
    return this->subdivideCellByCut(cf->getId(), plv0, plv1);
}

CellEdge *CellMesh::straightenCellEdgeSafe(CellEdge *ce)
{
    //Subdivide surrounding cells so that their borders will not be affected
    CellVertex *cv0 = (CellVertex *) ce->getVertex(0);
    CellVertex *cv1 = (CellVertex *) ce->getVertex(1);
    const int incidentFaceId0 = ce->getIncidentFaceId(0);
    const int incidentFaceId1 = ce->getIncidentFaceId(1);
    QVector<PolylineVertex *> subdivisionPLVs0;
    subdivisionPLVs0 << m_cellFaces[incidentFaceId0]->getPolyVertexCorrespondingToCellVertex(cv0)
                     << m_cellFaces[incidentFaceId0]->getPolyVertexCorrespondingToCellVertex(cv1);
    QVector<PolylineVertex *> subdivisionPLVs1; //could be PLVs on a neighboring cell edge (therefore not part of the other cell face)
    subdivisionPLVs1 << m_cellFaces[incidentFaceId1]->getPolyVertexCorrespondingToCellVertex(cv0)
                     << m_cellFaces[incidentFaceId1]->getPolyVertexCorrespondingToCellVertex(cv1);

    CellVertex *cvTemp0 = this->subdivideCellCentral(incidentFaceId0, subdivisionPLVs0);
    CellVertex *cvTemp1 = this->subdivideCellCentral(incidentFaceId1, subdivisionPLVs1);

    //Now the edge can safely be straightened
    CellEdge *newCE = this->straightenCellEdgeUnsafe(ce);

    //Merge surrounding cells again
    this->mergeCellsAroundVertex(cvTemp0);
    this->mergeCellsAroundVertex(cvTemp1);

    return newCE;
}

CellVertex *CellMesh::optimizeVertexWithSafetyCheck(CellVertex *cv)
{
    //check if boundary of surrounding macro cell has self intersections (parameterization cannot handle that)
    bool selfIntersection = false;
    QSet<Vertex *> cellVerticesOnBoundary;
    foreach (int incidentCfId, *cv->getIncidentFaceIds()) {
        CellFace *incidentCf = m_cellFaces[incidentCfId];
        const int numCV = incidentCf->getNumberOfVertices();
        for (int iCv = 0; iCv < numCV; iCv++) {
            CellVertex *incidentCv = (CellVertex *) incidentCf->getVertex(iCv);
            if (incidentCv != cv && !cellVerticesOnBoundary.contains(incidentCv))
                cellVerticesOnBoundary.insert(incidentCv);
            else if (cellVerticesOnBoundary.contains(incidentCv))
                selfIntersection = true;
        }
    }

    if (selfIntersection)
        return this->optimizeVertexSafe(cv);
    else
        return this->optimizeVertexUnsafe(cv);
}

CellVertex *CellMesh::optimizeVertexUnsafe(CellVertex *cv)
{
    QVector<CellVertex *> subdivideVertices;
    foreach (int eId, *cv->getIncidentEdgeIds())
        subdivideVertices.append((CellVertex *) m_cellEdges[eId]->getConnectedVertex(cv));

    CellFace *cfNew = this->mergeCellsAroundVertex(cv);

    QVector<PolylineVertex *> subdividePLVs;
    foreach (CellVertex *cvBorder, subdivideVertices)
        subdividePLVs.append(cfNew->getPolyVertexCorrespondingToCellVertex(cvBorder));

    return this->subdivideCellCentral(cfNew->getId(), subdividePLVs);
}

CellVertex *CellMesh::optimizeVertexSafe(CellVertex *cv)
{
    //subdivide surrounding edges
    QVector<CellVertex *> edgeMidpoints;
    const int numberOfIncidentCellEdges = cv->getIncidentEdgeIds()->size();
    for (int i = 0; i < numberOfIncidentCellEdges; i++) {
        CellEdge *ce = (CellEdge *) m_cellEdges[cv->getIncidentEdgeIds()->at(i)];
        //TODO it would be nice if the edge could get subdivided at some point t in [0,1]

        const int numPlv = ce->getPolylineVertices()->size();
        //TODO this check and fix could be moved to the cell edge class
        if (numPlv < 3) {//we need at least three PLVs to pick the center one
            //Each cell edge should have at least two PLV
            PolylineVertex *plv0 = ce->getPolylineVertices()->at(0);
            PolylineVertex *plv1 = ce->getPolylineVertices()->at(1);
            QVector3D helperPos = (plv0->getPosition() + plv1->getPosition())/2;
            PolylineVertex *helperPlv =  new PolylineVertex(helperPos, plv0->getSharedFaceId(plv1));

            ce->getPolylineVertices()->insert(1, helperPlv);

            //Reset parameterization of incident faces (need to be recomputed, since adding a vertex is not possible with the current implementation)
            m_cellFaceParameterizations[ce->getIncidentFaceId(0)] = 0;
            m_cellFaceParameterizations[ce->getIncidentFaceId(1)] = 0;
        }

        edgeMidpoints << this->subdivideEdge(ce, ce->getPolylineVertices()->at(numPlv/2));
    }

    //subdivide surrounding cells
    //QVector<CellVertex *> newCentralCVs;
    QVector<CellEdge *> temporaryCellEdges;
    for (int i = 0; i < numberOfIncidentCellEdges; i++) {
        const int cfId = cv->getIncidentFaceIds()->at(i);
        CellFace *cf = m_cellFaces[cfId];
        const int numIncidentVertices = cf->getNumberOfVertices();
        QVector<PolylineVertex *> subdivisionVertices;
        for (int i = 0; i < numIncidentVertices; i++) {
            if (edgeMidpoints.contains((CellVertex *) cf->getVertex(i)))
                    subdivisionVertices << cf->getPolyVertexCorrespondingToCellVertex((CellVertex *) cf->getVertex(i));
        }

        //newCentralCVs << this->subdivideCellCentral(cfId, subdivisionVertices);
        temporaryCellEdges << this->subdivideCellByCut(cfId, subdivisionVertices[0], subdivisionVertices[1]);
    }

    //Vertex is not safe to optimize
    CellVertex *optimizedCV = this->optimizeVertexUnsafe(cv);

    //Merge surrounding cells again
    for (int i = 0; i < numberOfIncidentCellEdges; i++)
        this->mergeCellsAroundEdge(temporaryCellEdges[i]);
        //this->mergeCellsAroundVertex(newCentralCVs[i]);

    //Merge subdivided edges again
    for (int i = 0; i < numberOfIncidentCellEdges; i++)
        this->mergeEdgesAroundUnecessaryVertex(edgeMidpoints[i]);

    //straighten affected edges
    foreach (int eId, *optimizedCV->getIncidentEdgeIds())
        this->straightenCellEdgeWithSafetyCheck(m_cellEdges[eId]);

    return optimizedCV;
}

void CellMesh::doSubdivisionForFullMesh()
{
    this->cleanUpMeshIndices(); //Necessary to avoid NULLs in vertex, edge, and face lists
    this->computeFullParameterization(m_parameterizationBorderType);    //Do this now, as this method can be run in parallel

    const int numberOfCellVertices = this->getNumberOfVertices();
    const int numberOfCellEdges = this->getNumberOfEdges();
    const int numberOfCellFaces = this->getNumberOfFaces();

    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = this->getEdge(i);
        this->subdivideEdge(ce, ce->getPolylineVertices()->at(ce->getPolylineVertices()->size()/2));
    }

    //TODO this can be run in parallel, but requires some additional work to be 100% thread safe
    for (int i = 0; i < numberOfCellFaces; i++) {
        QVector<PolylineVertex *> subdivideVertices;
        CellFace *cf = this->getFace(i);
        const int numVertices = cf->getNumberOfVertices();
        for (int j = 0; j < numVertices; j++)
            if (cf->getVertex(j)->getId() >= numberOfCellVertices) //that way, only newly generated vertices will be used (i.e. center points on edges)
                subdivideVertices.append(cf->getPolyVertexCorrespondingToCellVertex((CellVertex *) cf->getVertex(j)));

        this->subdivideCellCentral(i, subdivideVertices);
    }

    this->cleanUpMeshIndices();

    this->addInfo(QString("All cells subdivided. New number of cells: %1\n").arg(m_cellFaces.size()));
}

void CellMesh::doSubdivisionForOneCellFace(CellFace *cf)
{
    this->cleanUpMeshIndices();
    const int cfId = cf->getId();
    const int numberOfCellVertices = this->getNumberOfVertices();

    const int numberOfIncidentCellEdges = cf->getNumberOfEdges();
    QVector<CellEdge *> incidentEdges;  //Extract edges (cannot iterate over edges of cf, because subdivision adds new edges)
    for (int i = 0; i < numberOfIncidentCellEdges; i++)
        incidentEdges.append((CellEdge *) cf->getEdge(i));

    QVector<PolylineVertex *> dummy;
    for (int i = 0; i < numberOfIncidentCellEdges; i++) {
        CellEdge *ce = incidentEdges[i];
        this->subdivideEdge(ce, ce->getPolylineVertices()->at(ce->getPolylineVertices()->size()/2), &dummy);
    }

    QVector<PolylineVertex *> subdivideVertices;
    const int numVertices = cf->getNumberOfVertices();
    for (int j = 0; j < numVertices; j++)
        if (cf->getVertex(j)->getId() >= numberOfCellVertices)
            subdivideVertices.append(cf->getPolyVertexCorrespondingToCellVertex((CellVertex *) cf->getVertex(j)));

    this->subdivideCellCentral(cfId, subdivideVertices);

    this->cleanUpMeshIndices();

    this->addInfo(QString("Cell %1 subdivided. New number of cells: %2\n").arg(cfId).arg(m_cellFaces.size()));
}

void CellMesh::transformToDualMesh()
{
    this->cleanUpMeshIndices();

    const int numberOfInitialCellVertices = m_cellVertices.size();

    //Subdivide all cells (connect centers to centers of boundaries)
    this->doSubdivisionForFullMesh();

    //Merge cells around the vertices of the initial cellmesh
    for (int i = 0; i < numberOfInitialCellVertices; i++)
        this->mergeCellsAroundVertex(m_cellVertices[i]);

    this->cleanUpMeshIndices();

    //clean up degenerate cell vertices
    const int numberOfIntermediateCellVertices = m_cellVertices.size();
    for (int i = 0; i < numberOfIntermediateCellVertices; i++)
        if (m_cellVertices[i]->getIncidentEdgeIds()->size() == 2)
            this->mergeEdgesAroundUnecessaryVertex(m_cellVertices[i]);

    this->cleanUpMeshIndices();

    this->addInfo(QString("Dual mesh constructed. New number of cells: %1").arg(this->getNumberOfFaces()));
}

void CellMesh::computeFullParameterization(Parameterization::BorderType bt)
{
    this->cleanUpMeshIndices();
    const int numberOfCellFaces = m_cellFaces.size();

    #pragma omp parallel for
    for (int i = 0; i < numberOfCellFaces; i++)
        if (!m_cellFaceParameterizations[i] || m_cellFaceParameterizations[i]->getBorderType() != bt)
            m_cellFaceParameterizations[i] = Parameterization::harmonicMap(m_cellFaces[i], m_originalMesh, m_parameterizationWeightType, bt, m_parameterInterval, m_parameterizationStorageReservation_ExpectedEdgesPerVertex);
}

//void CellMesh::cutCellBy3DPlane(const int id, const QVector3D pointOnPlane, const QVector3D normal)
//{
//    CellFace *cf = m_cellFaces[id];

//    QVector<Vertex *> *meshVertices = cf->getMeshVertices();
//    QVector<Vertex *> meshVertices0;
//    QVector<Vertex *> meshVertices1;
//    QHash<Vertex *, int> vertexSideMap;

//    //First, subdivide the inner vertices
//    foreach (Vertex *v, *meshVertices) {
//        const double scalarProduct = QVector3D::dotProduct(v->getPosition() - pointOnPlane, normal);
//        if (scalarProduct > 0) {
//            meshVertices0.append(v);
//            vertexSideMap[v] = 0;
//        } else {
//            meshVertices1.append(v);
//            vertexSideMap[v] = 1;
//        }
//    }

//    //Also check on which side of the separating line the polyline vertices are located
//    const int numberOfCellEdges = cf->getNumberOfEdges();
//    for (int i = 0; i < numberOfCellEdges; i++) {
//        CellEdge *ce = (CellEdge *) cf->getEdge(i);
//        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
//            const double scalarProduct = QVector3D::dotProduct(plv->getPosition() - pointOnPlane, normal);
//            if (scalarProduct > 0)
//                vertexSideMap[plv] = 0;
//            else
//                vertexSideMap[plv] = 1;
//        }
//    }

//    //Second, subdivide the mesh edges. Add new polyline vertices on intersections with the separation line
//    QVector<Edge *> *meshEdges = cf->getMeshEdges();
//    QVector<Edge *> meshEdges0;
//    QVector<Edge *> meshEdges1;
//    QList<PolylineVertex *> polylineVertices;
//    const QVector3D normalizedNormal = normal.normalized();
//    foreach (Edge *e, *meshEdges) {
//        const int side0 = vertexSideMap[e->getVertex(0)];
//        const int side1 = vertexSideMap[e->getVertex(1)];
//        if (side0 == 0 && side1 == 0)
//            meshEdges0.append(e);
//        else if (side0 == 1 && side1 == 1)
//            meshEdges1.append(e);
//        else {
//            const QVector3D v0Pos = e->getVertex(0)->getPosition();
//            const QVector3D v1Pos = e->getVertex(1)->getPosition();
//            const double distToPlane = QVector3D::dotProduct(v0Pos - pointOnPlane, normalizedNormal);
//            QVector3D intersectionPos = v0Pos - distToPlane * normalizedNormal;
//            if (fabs(distToPlane) <= 0.01)
//                intersectionPos = v0Pos * 0.99 + v1Pos * 0.01;
//            PolylineVertex *plv = new PolylineVertex(intersectionPos, e);
//            polylineVertices.append(plv);
//        }
//    }

//    //Check for edges that are stored by the border and add polyline vertices if neccessary
//    QHash<Edge *, PolylineVertex *> crossingEdgePlvMap;
//    for (int i = 0; i < numberOfCellEdges; i++) {
//        CellEdge *ce = (CellEdge *) cf->getEdge(i);
//        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
//            if (plv->isMeshVertex()) {
//                foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
//                    Edge *e = m_originalMesh->getEdge(eId);
//                    if (!crossingEdgePlvMap.contains(e)) {
//                        crossingEdgePlvMap[e] = plv;

//                        Vertex *v = e->getConnectedVertex(plv->getMeshVertex());
//                        if (vertexSideMap.contains(v)) {
//                            const int side0 = vertexSideMap[plv];
//                            const int side1 = vertexSideMap[v];
//                            if (side0 != side1) {
//                                const QVector3D v0Pos = plv->getPosition();
//                                const QVector3D v1Pos = v->getPosition();
//                                const double distToPlane = QVector3D::dotProduct(v0Pos - pointOnPlane, normalizedNormal);
//                                QVector3D intersectionPos = v0Pos - distToPlane * normalizedNormal;
//                                if (fabs(distToPlane) <= 0.01)
//                                    intersectionPos = v0Pos * 0.99 + v1Pos * 0.01;

//                                PolylineVertex *plvNew = new PolylineVertex (intersectionPos, e);
//                                polylineVertices.append(plvNew);
//                            }
//                        }
////                    } else {
////                        PolylineVertex *otherPlv = crossingEdgePlvMap[e];
////                        const int side0 = vertexSideMap[plv];
////                        const int side1 = vertexSideMap[otherPlv];
////                        if (side0 != side1 && otherPlv->getPosition() != vSeparationStart->getPosition() && otherPlv->getPosition() != vSeparationEnd->getPosition()) {
////                            QVector2D intersectionParameter;
////                            PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, otherPlv, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
////                            polylineVertices.append(plvNew);
////                            newVertexParameters[plvNew] = intersectionParameter;
////                        }
//                    }
//                }
//            } else if (plv->isCrossingVertex() && plv->getPosition() != vSeparationStart->getPosition() && plv->getPosition() != vSeparationEnd->getPosition()) {
////                Edge *e = plv->getCrossingEdge();
////                if (!crossingEdgePlvMap.contains(e)) {
////                    crossingEdgePlvMap[e] = plv;

////                    for (int j = 0; j < 2; j++) {
////                        Vertex *v = e->getVertex(j);
////                        if (vertexSideMap.contains(v)) {
////                            const int side0 = vertexSideMap[plv];
////                            const int side1 = vertexSideMap[v];
////                            if (side0 != side1) {
////                                QVector2D intersectionParameter;
////                                PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, v, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
////                                polylineVertices.append(plvNew);
////                                newVertexParameters[plvNew] = intersectionParameter;
////                            }
////                        }
////                    }
////                } else {
////                    PolylineVertex *otherPlv = crossingEdgePlvMap[e];
////                    const int side0 = vertexSideMap[plv];
////                    const int side1 = vertexSideMap[otherPlv];
////                    if (side0 != side1 && otherPlv->getPosition() != vSeparationStart->getPosition() && otherPlv->getPosition() != vSeparationEnd->getPosition()) {
////                        QVector2D intersectionParameter;
////                        PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, otherPlv, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
////                        polylineVertices.append(plvNew);
////                        newVertexParameters[plvNew] = intersectionParameter;
////                    }
////                }
////            }
//        }
//    }


    //TODO:
    //-if cell edges exist: cut cell edges
    //-cut cell contents, extract polyline vertices for new cell edges
    //-sort polyline vertices to form a cell edge
    //-if cell edge has several connected components or more than two child cells would exist -> abort
//}

CellEdge *CellMesh::subdivideCellByCut(int id, PolylineVertex *vSeparationStart, PolylineVertex *vSeparationEnd)
{
    CellFace *cf = m_cellFaces[id];

    Parameterization *cellParameterization = this->getCellFaceParameterization(id);

    //Define cutting line
    const QVector2D projectedStart = cellParameterization->getParameter(vSeparationStart);
    const double projectedStartX = projectedStart.x();
    const double projectedStartY = projectedStart.y();
    const QVector2D projectedEnd = cellParameterization->getParameter(vSeparationEnd);
    const double projectedEndX = projectedEnd.x();
    const double projectedEndY = projectedEnd.y();

    const double dirX = projectedEndX - projectedStartX;
    const double dirY = projectedEndY - projectedStartY;

    const QVector3D normal = QVector3D::crossProduct(QVector3D(dirX, dirY, 0), QVector3D(0, 0, 1)).normalized();
    const QVector2D normal2D = normal.toVector2D();
    const double normalX = normal.x();
    const double normalY = normal.y();

    //Subdivide content of the cellface into two subsets
    QVector<Vertex *> *meshVertices = cf->getMeshVertices();
    QVector<Vertex *> meshVertices0;
    QVector<Vertex *> meshVertices1;
    QHash<Vertex *, int> vertexSideMap;

    //First, subdivide the inner vertices
    foreach (Vertex *v, *meshVertices) {
        const QVector2D paramValue = cellParameterization->getParameter(v);
        const double scalarProduct = ((paramValue.x()-projectedStartX)*normalX)+((paramValue.y()-projectedStartY)*normalY);
        if (scalarProduct > 0) {
            meshVertices0.append(v);
            vertexSideMap[v] = 0;
        } else {
            meshVertices1.append(v);
            vertexSideMap[v] = 1;
        }
    }

    //Also check on which side of the separating line the polyline vertices are located
    const int numberOfCellEdges = cf->getNumberOfEdges();
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            const QVector2D paramValue = cellParameterization->getParameter(plv);
            const double scalarProduct = ((paramValue.x()-projectedStartX)*normalX)+((paramValue.y()-projectedStartY)*normalY);
            if (scalarProduct > 0)
                vertexSideMap[plv] = 0;
            else
                vertexSideMap[plv] = 1;
        }
    }

    //Second, subdivide the mesh edges. Add new polyline vertices on intersections with the separation line
    QVector<Edge *> *meshEdges = cf->getMeshEdges();
    QVector<Edge *> meshEdges0;
    QVector<Edge *> meshEdges1;
    QList<PolylineVertex *> polylineVertices;
    QHash<PolylineVertex *, QVector2D> newVertexParameters;
    foreach (Edge *e, *meshEdges) {
        const int side0 = vertexSideMap[e->getVertex(0)];
        const int side1 = vertexSideMap[e->getVertex(1)];
        if (side0 == 0 && side1 == 0)
            meshEdges0.append(e);
        else if (side0 == 1 && side1 == 1)
            meshEdges1.append(e);
        else {
            QVector2D intersectionParameter;
            PolylineVertex *plv = new PolylineVertex(this->calculateIntersectionWithLine(e->getVertex(0), e->getVertex(1), cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
            polylineVertices.append(plv);
            newVertexParameters[plv] = intersectionParameter;
        }
    }

    //Check for edges that are stored by the border and add polyline vertices if neccessary
    QHash<Edge *, PolylineVertex *> crossingEdgePlvMap;
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            if (plv->isMeshVertex() && plv->getPosition() != vSeparationStart->getPosition() && plv->getPosition() != vSeparationEnd->getPosition()) {
                foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
                    Edge *e = m_originalMesh->getEdge(eId);
                    if (!crossingEdgePlvMap.contains(e)) {
                        crossingEdgePlvMap[e] = plv;

                        Vertex *v = e->getConnectedVertex(plv->getMeshVertex());
                        if (vertexSideMap.contains(v)) {
                            const int side0 = vertexSideMap[plv];
                            const int side1 = vertexSideMap[v];
                            if (side0 != side1) {
                                QVector2D intersectionParameter;
                                PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, v, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
                                polylineVertices.append(plvNew);
                                newVertexParameters[plvNew] = intersectionParameter;
                            }
                        }
                    } else {
                        PolylineVertex *otherPlv = crossingEdgePlvMap[e];
                        const int side0 = vertexSideMap[plv];
                        const int side1 = vertexSideMap[otherPlv];
                        if (side0 != side1 && otherPlv->getPosition() != vSeparationStart->getPosition() && otherPlv->getPosition() != vSeparationEnd->getPosition()) {
                            QVector2D intersectionParameter;
                            PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, otherPlv, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
                            polylineVertices.append(plvNew);
                            newVertexParameters[plvNew] = intersectionParameter;
                        }
                    }
                }
            } else if (plv->isCrossingVertex() && plv->getPosition() != vSeparationStart->getPosition() && plv->getPosition() != vSeparationEnd->getPosition()) {
                Edge *e = plv->getCrossingEdge();
                if (!crossingEdgePlvMap.contains(e)) {
                    crossingEdgePlvMap[e] = plv;

                    for (int j = 0; j < 2; j++) {
                        Vertex *v = e->getVertex(j);
                        if (vertexSideMap.contains(v)) {
                            const int side0 = vertexSideMap[plv];
                            const int side1 = vertexSideMap[v];
                            if (side0 != side1) {
                                QVector2D intersectionParameter;
                                PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, v, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
                                polylineVertices.append(plvNew);
                                newVertexParameters[plvNew] = intersectionParameter;
                            }
                        }
                    }
                } else {
                    PolylineVertex *otherPlv = crossingEdgePlvMap[e];
                    const int side0 = vertexSideMap[plv];
                    const int side1 = vertexSideMap[otherPlv];
                    if (side0 != side1 && otherPlv->getPosition() != vSeparationStart->getPosition() && otherPlv->getPosition() != vSeparationEnd->getPosition()) {
                        QVector2D intersectionParameter;
                        PolylineVertex *plvNew = new PolylineVertex (this->calculateIntersectionWithLine(plv, otherPlv, cellParameterization, projectedStart, normal2D, &intersectionParameter), e);
                        polylineVertices.append(plvNew);
                        newVertexParameters[plvNew] = intersectionParameter;
                    }
                }
            }
        }
    }

    //Sort list of new polyline vertices such that the first vertex is closest to vSeparationStart
    std::sort(polylineVertices.begin(), polylineVertices.end(), PolylineVertexComparator(&newVertexParameters, projectedStart));

    //Add start and end vertex to the list
    polylineVertices.prepend(new PolylineVertex(vSeparationStart));
    polylineVertices.append(new PolylineVertex(vSeparationEnd));

    //Check if input vertices are already cell vertices
    QVector<PolylineVertex *> vSplit(2);
    vSplit[0] = vSeparationStart;
    vSplit[1] = vSeparationEnd;
    CellVertex *cvSplit[2] = {0, 0};
    const int numberOfCellVertices = cf->getNumberOfVertices();
    for (int i = 0; i < numberOfCellVertices; i++) {
        if (cf->getVertex(i)->getPosition() == vSplit[0]->getPosition())
            cvSplit[0] = (CellVertex *) cf->getVertex(i);
        if (cf->getVertex(i)->getPosition() == vSplit[1]->getPosition())
            cvSplit[1] = (CellVertex *) cf->getVertex(i);
    }

    //If one of the vertices is no cell vertex -> split the containing edge at this vertex and create a new cell vertex
    for (int vId = 0; vId < 2; vId++) {
        if (!cvSplit[vId]) {
            PolylineVertex *splitVertex = vSplit[vId];

            //Find the cell edge that contains the split vertex
            CellEdge *edgeToSplit = 0;
            for (int i = 0; i < cf->getNumberOfEdges(); i++)
                if (((CellEdge *) cf->getEdge(i))->getPolylineVertices()->contains(splitVertex))
                    edgeToSplit = (CellEdge *) cf->getEdge(i);

            //Carefull: split vertex is destroyed after edge split -> vSplit will contain updated pointers
            cvSplit[vId] = this->subdivideEdge(edgeToSplit, splitVertex, &vSplit);
        }
    }   //At this point, edge splits are done and cell vertices are inserted as needed. We can continue as if no edge split was necessary
    //Assume now, that end vertices of the split edge are cvSplit[0 and 1]

    //Create new cell edge and update incident edge IDs of end vertices
    const int newCellFaceId = this->takeFreeFaceId();
    const int newCellEdgeId = this->takeFreeEdgeId();
    CellEdge *newCellEdge = new CellEdge(cvSplit[0], cvSplit[1], polylineVertices.toVector(), newCellEdgeId);
    newCellEdge->addIncidentFaceId(id);
    newCellEdge->addIncidentFaceId(newCellFaceId);
    m_cellEdges[newCellEdgeId] = newCellEdge;
    cvSplit[0]->addIncidentEdgeId(newCellEdgeId);
    cvSplit[1]->addIncidentEdgeId(newCellEdgeId);

    //Get cell edges of the current cell face (which should form a circle)
    QVector<CellEdge *> cEdges;
    for (int i = 0; i < cf->getNumberOfEdges(); i++)
        cEdges.append((CellEdge *) cf->getEdge(i));

    //Find the indices where the circle is split into two half-circles
    const int numberOfOriginalEdges = cEdges.size();
    int splitIndex0 = -1;
    int splitIndex1 = -1;
    CellVertex *sharedVertex = (CellVertex *) cEdges.last()->getSharedVertex(cEdges.first());
    if (cvSplit[0] == sharedVertex)
        splitIndex0 = 0;
    if (cvSplit[1] == sharedVertex)
        splitIndex1 = 0;
    for (int i = 1; i < numberOfOriginalEdges; i++) {
        sharedVertex = (CellVertex *) cEdges[i]->getSharedVertex(cEdges[i-1]);
        if (cvSplit[0] == sharedVertex)
            splitIndex0 = i;
        if (cvSplit[1] == sharedVertex)
            splitIndex1 = i;
    }

    if (splitIndex0 > splitIndex1) {
        int tmp = splitIndex0;
        splitIndex0 = splitIndex1;
        splitIndex1 = tmp;
    }

    //Subdivide circle of edges into two half circles
    QVector<Edge *> cEdgesSubdivided[2] = {QVector<Edge *>(), QVector<Edge *>()};
    for (int i = splitIndex0; i < splitIndex1; i++)
        cEdgesSubdivided[0].append(cEdges[i]);
    for (int i = splitIndex1; i < splitIndex0 + numberOfOriginalEdges; i++)
        cEdgesSubdivided[1].append(cEdges[i % numberOfOriginalEdges]);

    //Close half circles
    cEdgesSubdivided[0].append(newCellEdge);
    cEdgesSubdivided[1].append(newCellEdge);

    //Extract vertices of new circles
    QVector<Vertex *> cVerticesSubdivided[2] = {QVector<Vertex *>(), QVector<Vertex *>()};
    for (int i = 0; i < 2; i++) {
        cVerticesSubdivided[i].append(cEdgesSubdivided[i].last()->getSharedVertex(cEdgesSubdivided[i].first()));
        for (int j = 1; j < cEdgesSubdivided[i].size(); j++)
            cVerticesSubdivided[i].append(cEdgesSubdivided[i].at(j-1)->getSharedVertex(cEdgesSubdivided[i].at(j)));
        if (cVerticesSubdivided[i].size() < 3) {
            qDebug() <<  "Error in subdivision by cut: Not enough cell vertices in new cell (< 3)";
            exit(-1);
        }
    }

    //Match new cell edge circles to the contained mesh faces
    //Take some a cell vertex that is not one of the newly created cell vertices
    Vertex *cvCompare = cVerticesSubdivided[0].at(1);
    QVector2D paramValue = cellParameterization->getParameter(cvCompare);
    double scalarProduct = ((paramValue.x()-projectedStartX)*normalX)+((paramValue.y()-projectedStartY)*normalY);

    //Check on which side of the separating line the vertex is located
    int cf0Edges = -1;
    if (scalarProduct > 0)
        cf0Edges = 0;
    else
        cf0Edges = 1;

    //Create new cell faces
    CellFace *cf0 = new CellFace(cEdgesSubdivided[cf0Edges], meshEdges0, meshVertices0, newCellFaceId);
    CellFace *cf1 = new CellFace(cEdgesSubdivided[1 - cf0Edges], meshEdges1, meshVertices1, id);
    m_cellFaces[newCellFaceId] = cf0;
    m_cellFaces[id] = cf1;
    m_cellFaceParameterizations[id] = 0;
    m_cellFaceParameterizations[newCellFaceId] = 0;

    //Update incident face IDs of cell edges (only for cf0)
    for (int i = 0; i < cEdgesSubdivided[cf0Edges].size() - 1; i++)     //skip the last cell edge, since it is the newly created edge
        ((CellEdge *) cEdgesSubdivided[cf0Edges].at(i))->updateIncidentFaceId(id, newCellFaceId);

    //Update incident face IDs of cell vertices (only for cf0)
    cvSplit[0]->addIncidentFaceId(newCellFaceId);
    cvSplit[1]->addIncidentFaceId(newCellFaceId);
    for (int i = 1; i < cVerticesSubdivided[cf0Edges].size()-1; i++)    //skip the first and last cell vertex (endpoints of new cell edge)
        ((CellVertex *) cVerticesSubdivided[cf0Edges].at(i))->updateIncidentFaceId(id, newCellFaceId);

    delete cellParameterization;
    delete cf;

    return newCellEdge;
}

CellVertex *CellMesh::subdivideCellCentral(int id, QVector<PolylineVertex *> subdivisionVertices)
{
    CellFace *cf = m_cellFaces[id];
    const int numberOfChildCells = subdivisionVertices.size();
    const int numberOfCellEdges = cf->getNumberOfEdges();
//    if (numberOfChildCells < 3) {
//        qDebug() << "Error in central subdivision: At least three separation vertices are required! Cell:" << id;
//        exit(-2);
//    }

    //Compute parameterization
    Parameterization *cellParameterization = this->getCellFaceParameterization(id);

    //Parameterize separation vertices
    QVector<double> separationAnglesUnsorted(numberOfChildCells);
    QVector<QVector2D> separationCoordinatesUnsorted(numberOfChildCells);
    for (int i = 0; i < numberOfChildCells; i++) {
        separationAnglesUnsorted[i] = cellParameterization->getPolarAngle(subdivisionVertices[i]);
        separationCoordinatesUnsorted[i] = cellParameterization->getParameter(subdivisionVertices[i]);
    }

    //Sort subdivision vertices
    QVector<PolylineVertex *> tmpList(subdivisionVertices);
    QVector<double> separationAngles;
    QVector<QVector2D> separationCoordinates;
    QVector<PolylineVertex *> separationVertices;
    while (!tmpList.isEmpty()) {
        int index = 0;
        double lowest = separationAnglesUnsorted[0];
        for (int i = 1; i < separationAnglesUnsorted.size(); i++) {
            if (separationAnglesUnsorted[i] < lowest) {
                index = i;
                lowest = separationAnglesUnsorted[i];
            }
        }
        separationAngles.append(separationAnglesUnsorted[index]);
        separationCoordinates.append(separationCoordinatesUnsorted[index]);
        separationVertices.append(tmpList[index]);

        separationAnglesUnsorted.remove(index);
        separationCoordinatesUnsorted.remove(index);
        tmpList.remove(index);
    }

    //subdivide inner vertices
    QVector<Vertex *> *origMeshVertices = cf->getMeshVertices();
    QVector<QVector<Vertex *> > meshVertices(numberOfChildCells);
    QHash<Vertex *, int> vertexSliceMap;

    //Child cell i contains the vertices between separation vertices i and i+1
    foreach (Vertex *v, *origMeshVertices) {
        const int index = this->findSlice(cellParameterization->getPolarAngle(v), separationAngles);
        meshVertices[index].append(v);
        vertexSliceMap[v] = index;
    }

    //Also check on which side of the separating lines the polyline vertices are located
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            const int index = this->findSlice(cellParameterization->getPolarAngle(plv), separationAngles);
            vertexSliceMap[plv] = index;
        }
    }

    //Now, subdivide the mesh edges. Add new polyline vertices on intersections with the separation lines
    QVector<Edge *> *origMeshEdges = cf->getMeshEdges();
    QVector<QVector<Edge *> > meshEdges(numberOfChildCells);
    QVector<QList<PolylineVertex *> > polylineVertices(numberOfChildCells);
    QHash<PolylineVertex *, QVector2D> newVertexParameters;
    const double centerParameterValue = m_parameterInterval.first + (m_parameterInterval.second - m_parameterInterval.first)/2;
    const QVector2D origin(centerParameterValue, centerParameterValue);     //in theory, only this value has to be updated if we want to do subdivision on one of the cell edges

    const float minLambda0 = 0.00001;
    const float maxLambda0 = 0.99999;
    const float minLambda1 = -1;
    const float maxLambda1 = 0.99999;

    foreach (Edge *e, *origMeshEdges) {
        const int slice0 = vertexSliceMap[e->getVertex(0)];
        const int slice1 = vertexSliceMap[e->getVertex(1)];
        if (slice0 == slice1) {
            meshEdges[slice0].append(e);
        } else {
            QVector2D param0 = cellParameterization->getParameter(e->getVertex(0));
            QVector2D param1 = cellParameterization->getParameter(e->getVertex(1));
            for (int i = 0; i < numberOfChildCells; i++) {
                QVector2D lambda = this->calculateIntersectionCoefficients(param0, param1, origin, separationCoordinates[i]);

                if (lambda.y() <= 1 && lambda.x() >= 0 && lambda.x() <= 1) {
                    //Ensure that the intersection has a minimal distance to the endpoints
                    this->restrictToInterval(&lambda, minLambda0, maxLambda0, minLambda1, maxLambda1);

                    //QVector2D intersectionParam = (1 - lambda.y()) * separationCoordinates[i];   // + lambda * (0, 0)
                    QVector2D intersectionParam = lambda.x() * param0 + (1 - lambda.x()) * param1;
                    //... or intersectionParam = lambda.x() * param0 + (1 - lambda.x()) * param1
                    QVector3D intersection3D = lambda.x() * e->getVertex(0)->getPosition() + (1 - lambda.x()) * e->getVertex(1)->getPosition();

                    PolylineVertex *plv = new PolylineVertex(intersection3D, e);
                    polylineVertices[i].append(plv);
                    newVertexParameters[plv] = intersectionParam;
                }
            }
        }
    }

    //Check for edges that are stored by the border and add polyline vertices if neccessary
    QHash<Edge *, PolylineVertex *> crossingEdgePlvMap;
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            if (plv->isMeshVertex()) {
                foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
                    Edge *e = m_originalMesh->getEdge(eId);
                    if (!crossingEdgePlvMap.contains(e)) {
                        crossingEdgePlvMap[e] = plv;

                        Vertex *v = e->getConnectedVertex(plv->getMeshVertex());
                        if (vertexSliceMap.contains(v)) {
                            const int slice0 = vertexSliceMap[plv];
                            const int slice1 = vertexSliceMap[v];
                            if (slice0 != slice1) {
                                QVector2D param0 = cellParameterization->getParameter(plv);
                                QVector2D param1 = cellParameterization->getParameter(v);
                                for (int i = 0; i < numberOfChildCells; i++) {
                                    QVector2D lambda = this->calculateIntersectionCoefficients(param0, param1, origin, separationCoordinates[i]);
                                    if (lambda.y() <= 1 && lambda.x() >= minLambda0 && lambda.x() <= maxLambda0) {
                                        this->restrictToInterval(&lambda, minLambda0, maxLambda0, minLambda1, maxLambda1);
                                        QVector2D intersectionParam = (1 - lambda.y()) * separationCoordinates[i];
                                        QVector3D intersection3D = lambda.x() * plv->getPosition() + (1 - lambda.x()) * v->getPosition();
                                        PolylineVertex *plvNew = new PolylineVertex(intersection3D, e);
                                        polylineVertices[i].append(plvNew);
                                        newVertexParameters[plvNew] = intersectionParam;
                                    }
                                }
                            }
                        }
                    } else {
                        PolylineVertex *otherPlv = crossingEdgePlvMap[e];
                        const int slice0 = vertexSliceMap[plv];
                        const int slice1 = vertexSliceMap[otherPlv];

                        if (slice0 != slice1) {
                            QVector2D param0 = cellParameterization->getParameter(plv);
                            QVector2D param1 = cellParameterization->getParameter(otherPlv);
                            for (int i = 0; i < numberOfChildCells; i++) {
                                QVector2D lambda = this->calculateIntersectionCoefficients(param0, param1, origin, separationCoordinates[i]);
                                if (lambda.y() <= 1 && lambda.x() >= minLambda0 && lambda.x() <= maxLambda0) {
                                    this->restrictToInterval(&lambda, minLambda0, maxLambda0, minLambda1, maxLambda1);
                                    QVector2D intersectionParam = (1 - lambda.y()) * separationCoordinates[i];
                                    QVector3D intersection3D = lambda.x() * plv->getPosition() + (1 - lambda.x()) * otherPlv->getPosition();
                                    PolylineVertex *plvNew = new PolylineVertex(intersection3D, e);
                                    polylineVertices[i].append(plvNew);
                                    newVertexParameters[plvNew] = intersectionParam;
                                }
                            }
                        }
                    }
                }
            } else if (plv->isCrossingVertex()) {
                Edge *e = plv->getCrossingEdge();
                if (!crossingEdgePlvMap.contains(e)) {
                    crossingEdgePlvMap[e] = plv;

                    for (int j = 0; j < 2; j++) {
                        Vertex *v = e->getVertex(j);
                        if (vertexSliceMap.contains(v)) {
                            const int slice0 = vertexSliceMap[plv];
                            const int slice1 = vertexSliceMap[v];
                            if (slice0 != slice1) {
                                QVector2D param0 = cellParameterization->getParameter(plv);
                                QVector2D param1 = cellParameterization->getParameter(v);
                                for (int i = 0; i < numberOfChildCells; i++) {
                                    QVector2D lambda = this->calculateIntersectionCoefficients(param0, param1, origin, separationCoordinates[i]);
                                    if (lambda.y() <= 1 && lambda.x() >= minLambda0 && lambda.x() <= maxLambda0) {
                                        this->restrictToInterval(&lambda, minLambda0, maxLambda0, minLambda1, maxLambda1);
                                        QVector2D intersectionParam = (1 - lambda.y()) * separationCoordinates[i];
                                        QVector3D intersection3D = lambda.x() * plv->getPosition() + (1 - lambda.x()) * v->getPosition();
                                        PolylineVertex *plvNew = new PolylineVertex(intersection3D, e);
                                        polylineVertices[i].append(plvNew);
                                        newVertexParameters[plvNew] = intersectionParam;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    PolylineVertex *otherPlv = crossingEdgePlvMap[e];
                    const int slice0 = vertexSliceMap[plv];
                    const int slice1 = vertexSliceMap[otherPlv];
                    if (slice0 != slice1) {
                        QVector2D param0 = cellParameterization->getParameter(plv);
                        QVector2D param1 = cellParameterization->getParameter(otherPlv);
                        for (int i = 0; i < numberOfChildCells; i++) {
                            QVector2D lambda = this->calculateIntersectionCoefficients(param0, param1, origin, separationCoordinates[i]);
                            if (lambda.y() <= 1 && lambda.x() >= minLambda0 && lambda.x() <= maxLambda0) {
                                this->restrictToInterval(&lambda, minLambda0, maxLambda0, minLambda1, maxLambda1);
                                QVector2D intersectionParam = (1 - lambda.y()) * separationCoordinates[i];
                                QVector3D intersection3D = lambda.x() * plv->getPosition() + (1 - lambda.x()) * otherPlv->getPosition();
                                PolylineVertex *plvNew = new PolylineVertex(intersection3D, e);
                                polylineVertices[i].append(plvNew);
                                newVertexParameters[plvNew] = intersectionParam;
                            }
                        }
                    }
                }
            }
        }
    }

    //Sort lists of new polyline vertices such that the first vertex is closest to (0, 0)
    for (int i = 0; i < numberOfChildCells; i++) {
        std::sort(polylineVertices[i].begin(), polylineVertices[i].end(), PolylineVertexComparator(&newVertexParameters, origin));
        polylineVertices[i].append(new PolylineVertex(separationVertices[i]));
    }

    //Find face that contains the origin
    QVector<int> incidentMeshFaces;
    incidentMeshFaces.append(polylineVertices.first().first()->getCrossingEdge()->getIncidentFaceId(0));
    incidentMeshFaces.append(polylineVertices.first().first()->getCrossingEdge()->getIncidentFaceId(1));
    for (int i = 1; i < numberOfChildCells; i++) {
        int fId0 = polylineVertices.at(i).at(0)->getCrossingEdge()->getIncidentFaceId(0);
        int fId1 = polylineVertices.at(i).at(0)->getCrossingEdge()->getIncidentFaceId(1);
        for (int j = incidentMeshFaces.size() - 1; j >= 0 ; j--) {  //Go backwards to avoid problems with indices after removing
            if (incidentMeshFaces[j] != fId0 && incidentMeshFaces[j] != fId1)
                incidentMeshFaces.remove(j);
        }
    }

    //Determine the center vertex
    PolylineVertex *plvCenter = 0;
    if (!incidentMeshFaces.isEmpty()) {
        //Take first three vertices of face and calculate parameters for for linear combination of origin
        Face *centerFace = m_originalMesh->getFace(incidentMeshFaces.first());
        Vertex *centerVertices[3] = {centerFace->getVertex(0),
                                     centerFace->getVertex(1),
                                     centerFace->getVertex(2)};
        QVector2D paramValues[3] = {cellParameterization->getParameter(centerVertices[0]),
                                    cellParameterization->getParameter(centerVertices[1]),
                                    cellParameterization->getParameter(centerVertices[2])};

        //barycentric coordinates (wikipedia)
        //const QVector2D vec02 = paramValues[0] - paramValues[2];
        const QVector2D vec02 = paramValues[0] - paramValues[2];
        const QVector2D vec12 = paramValues[1] - paramValues[2];
        const double det = (vec02.x() * vec12.y()) - (vec02.y() * vec12.x());
        if (det != 0) {
            const double l0 = (vec12.y() * (centerParameterValue - paramValues[2].x()) - vec12.x() * (centerParameterValue - paramValues[2].y()))/det;
            const double l1 = (-vec02.y() * (centerParameterValue - paramValues[2].x()) + vec02.x() * (centerParameterValue - paramValues[2].y()))/det;
            const double l2 = 1 - l1 - l0;
            const QVector3D posCenter = l0 * centerVertices[0]->getPosition() + l1 * centerVertices[1]->getPosition() + l2 * centerVertices[2]->getPosition();
            plvCenter = new PolylineVertex(posCenter, centerFace->getId());
        }
        //Take center of incident face as center vertex
//        QVector3D posCenter(0, 0, 0);
//        Face *centerFace = m_originalMesh->getFace(incidentMeshFaces.first());
//        for (int i = 0; i < centerFace->getNumberOfVertices(); i++)
//            posCenter += centerFace->getVertex(i)->getPosition();

//        posCenter /= centerFace->getNumberOfVertices();
//        plvCenter = new PolylineVertex(posCenter);
    }

    if (incidentMeshFaces.isEmpty() || !plvCenter) {
        qDebug() << "Problem in central subdivision: Inner polyline vertices do not seem to share a face or shared face does not contain origin";
        //Alternative approach: Find polyline Vertex that is closest to the center and add it as new center
        plvCenter = polylineVertices[0][0];
        double distToCenter = newVertexParameters[plvCenter].length();
        for (int i = 1; i < numberOfChildCells; i++) {
            PolylineVertex *plv = polylineVertices[i][0];
            double newDist = newVertexParameters[plv].length();
            if (newDist < distToCenter) {
                plvCenter = plv;
                distToCenter = newDist;
            }
        }
    }

    //Add copy of central vertex to each list of polyline vertices
    for (int i = 0; i < numberOfChildCells; i++)
        if (polylineVertices[i][0] != plvCenter)
            polylineVertices[i].prepend(new PolylineVertex(plvCenter));

    //Get new face IDs for the child cells
    QVector<int> newCellFaceIds(numberOfChildCells, -1);
    for (int i = 0; i < numberOfChildCells - 1; i++)
        newCellFaceIds[i] = this->takeFreeFaceId();
    newCellFaceIds.last() = id;

    //Create new Cell Vertex for central subdivision
    int newVertexId = this->takeFreeVertexId();
    CellVertex *cvCentral = new CellVertex(plvCenter->getPosition(), newVertexId);
    m_cellVertices[newVertexId] = cvCentral;
    for (int i = 0; i < numberOfChildCells; i++)
        cvCentral->addIncidentFaceId(newCellFaceIds[i]);


    //Check if input vertices are already cell vertices
    QVector<CellVertex *> cvSplit(numberOfChildCells, 0);
    const int numberOfCellVertices = cf->getNumberOfVertices();
    for (int i = 0; i < numberOfCellVertices; i++)
        for (int j = 0; j < numberOfChildCells; j++)
            if (cf->getVertex(i)->getPosition() == separationVertices[j]->getPosition())
                cvSplit[j] = (CellVertex *) cf->getVertex(i);

    //If one of the vertices is no cell vertex -> split the containing edge at this vertex and create a new cell vertex
    QVector<PolylineVertex *> vSplit(separationVertices);
    for (int vId = 0; vId < numberOfChildCells; vId++) {
        if (!cvSplit[vId]) {
            PolylineVertex *splitVertex = vSplit[vId];

            //Find the cell edge that contains the split vertex
            CellEdge *edgeToSplit = 0;
            for (int i = 0; i < cf->getNumberOfEdges(); i++)
                if (((CellEdge *) cf->getEdge(i))->getPolylineVertices()->contains(splitVertex))
                    edgeToSplit = (CellEdge *) cf->getEdge(i);

            //Carefull: split vertex is destroyed after edge split -> vSplit will contain updated pointers
            cvSplit[vId] = this->subdivideEdge(edgeToSplit, splitVertex, &vSplit);
        }
    }   //At this point, edge splits are done and cell vertices are inserted as needed. We can continue as if no edge split was necessary

    //Create new cell edges and update incident edge IDs of end vertices
    QVector<CellEdge *> innerCellEdges(numberOfChildCells);
    for (int i = 0; i < numberOfChildCells; i++) {
        const int newCellEdgeId = this->takeFreeEdgeId();
        CellEdge *newCellEdge = new CellEdge(cvCentral, cvSplit[i], polylineVertices[i].toVector(), newCellEdgeId);
        newCellEdge->addIncidentFaceId(newCellFaceIds[i]);
        newCellEdge->addIncidentFaceId(newCellFaceIds[(i + numberOfChildCells - 1) % numberOfChildCells]);
        m_cellEdges[newCellEdgeId] = newCellEdge;
        cvCentral->addIncidentEdgeId(newCellEdgeId);
        cvSplit[i]->addIncidentEdgeId(newCellEdgeId);

        innerCellEdges[i] = newCellEdge;
    }

    //Find the cell edges and vertices that contain each child cell
    QVector<QVector<Edge *> > cEdges(numberOfChildCells);
    QVector<QVector<CellVertex *> > outerCellVertices(numberOfChildCells);
    const int numberOfCellEdgesAfterSubdivision = cf->getNumberOfEdges();
    for (int i = 0; i < numberOfChildCells; i++) {
        CellEdge *ceSlice0 = innerCellEdges[i];
        CellEdge *ceSlice1 = innerCellEdges[(i + 1) % numberOfChildCells];
        cEdges[i].append(ceSlice1);
        cEdges[i].append(ceSlice0);

        //Find Cell vertex where the outer edges start -> index of first edge corresponds to index of vertex
        CellVertex *cv0 = (CellVertex *) ceSlice0->getConnectedVertex(cvCentral);
        CellVertex *cv1 = (CellVertex *) ceSlice1->getConnectedVertex(cvCentral);

        int index = 0;
        while (cf->getVertex(index) != cv0)
            index++;

        while (cf->getVertex(index) != cv1) {
            outerCellVertices[i].append((CellVertex *) cf->getVertex(index));
            cEdges[i].append(cf->getEdge(index));
            index = (index + 1) % numberOfCellEdgesAfterSubdivision;
        }
        outerCellVertices[i].append(cv1);
    }

    //Create new cell faces
    for (int i = 0; i < numberOfChildCells; i++) {
        int newId = newCellFaceIds[i];
        CellFace *cfNew = new CellFace(cEdges[i], meshEdges[i], meshVertices[i], newId);
        m_cellFaces[newId] = cfNew;
        m_cellFaceParameterizations[newId] = 0;
        const int numEdges = cEdges[i].size();

        //Update incident edge IDs of cell edges
        for (int j = 2; j < numEdges; j++)      //Only update the outer edges (inner edges are already done)
            ((CellEdge *) cEdges[i][j])->updateIncidentFaceId(id, newId);

        //Update incident face IDs of cell vertices (except center vertex)
        int numOuterVertices = outerCellVertices[i].size();
        for (int j = 0; j < numOuterVertices - 1; j++)
            outerCellVertices[i][j]->updateIncidentFaceId(id, newId);
        outerCellVertices.at(i).last()->addIncidentFaceId(newId);   //Just add the new ID (update will be done in next loop)
    }

    delete cellParameterization;
    delete cf;

    return cvCentral;

    //For debuging
//    for (int i = 0; i < numberOfChildCells; i++) {
//        newVertexParameters[polylineVertices.at(i).first()] = QVector2D(0, 0);
//        newVertexParameters[polylineVertices.at(i).last()] = cellParameterization->getParameter(separationVertices[i]);
//    }
//    cellParameterization->saveAsImageWithSubdivisionLines(polylineVertices, newVertexParameters, QString("../../img_%1.png").arg(id));
}

CellFace *CellMesh::mergeCellsAroundEdge(CellEdge *separatingEdge)
{
    //Get incident faces
    int cfId0 = separatingEdge->getIncidentFaceId(0);
    int cfId1 = separatingEdge->getIncidentFaceId(1);
    if (cfId0 > cfId1) {
        int tmp = cfId0;
        cfId0 = cfId1;
        cfId1 = tmp;
    }

    CellFace *cf0 = m_cellFaces[cfId0];
    CellFace *cf1 = m_cellFaces[cfId1];

    //Get the cell edges that form a circle around the union of cf0 and cf1 (i.e. remove the "merge" edge and glue together the two half circles)
    const int numberOfEdges0 = cf0->getNumberOfEdges();
    const int numberOfEdges1 = cf1->getNumberOfEdges();

    int cf0EdgeIndex = 0;
    while (cf0->getEdge(cf0EdgeIndex) != separatingEdge)
        cf0EdgeIndex++;

    int cf1EdgeIndex = 0;
    while (cf1->getEdge(cf1EdgeIndex) != separatingEdge)
        cf1EdgeIndex++;

    QVector<Edge *> halfCircle0;
    for (int i = cf0EdgeIndex + 1; i < numberOfEdges0; i++)
        halfCircle0.append(cf0->getEdge(i));
    for (int i = 0; i < cf0EdgeIndex; i++)
        halfCircle0.append(cf0->getEdge(i));

    QVector<Edge *> halfCircle1;
    for (int i = cf1EdgeIndex + 1; i < numberOfEdges1; i++)
        halfCircle1.append(cf1->getEdge(i));
    for (int i = 0; i < cf1EdgeIndex; i++)
        halfCircle1.append(cf1->getEdge(i));

    QVector<Edge *> totalCircle(halfCircle0);
    if (halfCircle0.last()->getSharedVertex(halfCircle1.first()))
        totalCircle << halfCircle1;
    else
        for (int i = halfCircle1.size() - 1; i >= 0; i--)
            totalCircle.append(halfCircle1[i]);

    //Get vertices and edges contained in the two faces
    QVector<Vertex *> meshVertices(*cf0->getMeshVertices());
    meshVertices << *cf1->getMeshVertices();

    QVector<Edge *> meshEdges(*cf0->getMeshEdges());
    meshEdges << *cf1->getMeshEdges();

    //Get vertices and edges stored by the cell edge (skip the endpoints -> their information is redundant to the endpoints of incident cell edges)
    const int numberOfPolyVertices = separatingEdge->getPolylineVertices()->size();
    for (int i = 1; i < numberOfPolyVertices - 1; i++) {
        PolylineVertex *plv = separatingEdge->getPolylineVertices()->at(i);
        if (plv->isMeshVertex())
            meshVertices.append(plv->getMeshVertex());
    }

    //Second parse is neccessary (we first need to extract all mesh vertices, before checking the mesh edges)
    for (int i = 1; i < numberOfPolyVertices - 1; i++) {
        PolylineVertex *plv = separatingEdge->getPolylineVertices()->at(i);
        if (plv->isMeshVertex()) {
            Vertex *v = plv->getMeshVertex();
            foreach (int eId, *v->getIncidentEdgeIds()) {
                Edge *e = m_originalMesh->getEdge(eId);
                if (meshVertices.contains(e->getConnectedVertex(v)) && !meshEdges.contains(e))
                    meshEdges.append(e);
            }
        } else if (plv->isCrossingVertex()) {
            Edge *e = plv->getCrossingEdge();
            if (meshVertices.contains(e->getVertex(0)) && meshVertices.contains(e->getVertex(1)) && !meshEdges.contains(e))
                meshEdges.append(e);
        }
    }

    //Create new cell face
    CellFace *cfNew = new CellFace(totalCircle, meshEdges, meshVertices, cfId0);
    m_cellFaces[cfId0] = cfNew;

    //Update incident face ids of adjacent edges
    for (int i = 0; i < numberOfEdges1-1; i++)
        ((CellEdge *) halfCircle1[i])->updateIncidentFaceId(cfId1, cfId0);

    //Update incident face ids of adjacent vertices
    for (int i = 1; i < numberOfEdges1-1; i++)
        ((CellVertex *) halfCircle1[i]->getSharedVertex(halfCircle1[i-1]))->updateIncidentFaceId(cfId1, cfId0);

    //Update incident edge and face ids of the end vertices
    const int sepEdgeId = separatingEdge->getId();
    for (int i = 0; i < 2; i++) {
        CellVertex *cv = (CellVertex *) separatingEdge->getVertex(i);
        cv->updateIncidentEdgeId(sepEdgeId, -1);
        cv->updateIncidentFaceId(cfId1, -1);
    }

    //Delete old faces and edge and free up their indices
    m_freeFaceIds.append(cfId1);
    m_freeEdgeIds.append(sepEdgeId);
    if (m_cellFaceParameterizations[cfId0])
        delete m_cellFaceParameterizations[cfId0];
    if (m_cellFaceParameterizations[cfId1])
        delete m_cellFaceParameterizations[cfId1];
    m_cellFaceParameterizations[cfId0] = 0;
    m_cellFaceParameterizations[cfId1] = 0;

    m_cellFaces[cfId1] = 0;
    m_cellEdges[separatingEdge->getId()] = 0;
    delete cf0;
    delete cf1;
    delete separatingEdge;

    return cfNew;
}

CellFace *CellMesh::mergeCellsAroundVertex(CellVertex *centerVertex)
{
    //Extract incident faces and edges
    QVector<int> *incidentFaceIds = centerVertex->getIncidentFaceIds();
    const int numberOfIncidentFaces = incidentFaceIds->size();
    QVector<CellFace *> incidentFaces(numberOfIncidentFaces);
    for (int i = 0; i < numberOfIncidentFaces; i++)
        incidentFaces[i] = m_cellFaces[incidentFaceIds->at(i)];

    const int numberOfIncidentEdges = centerVertex->getIncidentEdgeIds()->size();
    QVector<CellEdge *> incidentEdges(numberOfIncidentEdges);
    for (int i = 0; i < numberOfIncidentEdges; i++)
        incidentEdges[i] = m_cellEdges[centerVertex->getIncidentEdgeIds()->at(i)];

    //Find the cell edges that form a circle around the union of the incident cell faces
    QVector<Edge *> borderCellEdges;
    CellEdge *ceFirst = m_cellEdges[centerVertex->getIncidentEdgeIds()->at(0)];
    int nextFaceId = ceFirst->getIncidentFaceId(0);
    CellEdge *ceStart = ceFirst;
    CellEdge *ceNext = 0;

    QVector<int> faceIdsForUpdate;
    while (ceNext != ceFirst) {
        CellFace *cf = m_cellFaces[nextFaceId];
        const int numEdgesInFace = cf->getNumberOfEdges();
        int eId = cf->getLocalEdgeId(ceStart);

        int dir = 1;
        if (cf->getEdge((eId + dir) % numEdgesInFace)->getSharedVertex(ceStart) == centerVertex)
            dir *= -1;
        CellEdge *ceEnd = (CellEdge *) cf->getEdge((eId + numEdgesInFace - dir) % numEdgesInFace);
        eId = (eId + numEdgesInFace + dir) % numEdgesInFace;
        ceNext = (CellEdge *) cf->getEdge(eId);
        while (ceNext != ceEnd) {
            borderCellEdges.append(ceNext);
            faceIdsForUpdate.append(cf->getId());
            eId = (eId + numEdgesInFace + dir) % numEdgesInFace;
            ceNext = (CellEdge *) cf->getEdge(eId);
        }
        nextFaceId = ceNext->getIncidentFaceId(0);
        if (nextFaceId == cf->getId())
            nextFaceId = ceNext->getIncidentFaceId(1);
        ceStart = ceNext;
    }

    //Get vertices and edges contained in the incident faces
    QVector<Vertex *> meshVertices;
    QVector<Edge *> meshEdges;
    for (int i = 0; i < numberOfIncidentFaces; i++) {
        meshVertices << *incidentFaces.at(i)->getMeshVertices();
        meshEdges << *incidentFaces.at(i)->getMeshEdges();
    }

    //Get an arbitrary copy of the center vertex and extract its information
    PolylineVertex *plvCenter = incidentFaces.first()->getPolyVertexCorrespondingToCellVertex(centerVertex);
    if (plvCenter->isMeshVertex())
        meshVertices.append(plvCenter->getMeshVertex());

    //Get vertices and edges stored by the inner cell edges (skip the endpoints -> their information is redundant to the endpoints of incident border cell edges)
    foreach (CellEdge *ce, incidentEdges) {
        const int numberOfPolyVertices = ce->getPolylineVertices()->size();
        for (int i = 1; i < numberOfPolyVertices - 1; i++) {
            PolylineVertex *plv = ce->getPolylineVertices()->at(i);
            if (plv->isMeshVertex())
                meshVertices.append(plv->getMeshVertex());
        }
    }

    //Second parse is neccessary (we first need to extract all mesh vertices, before checking the mesh edges)
    if (plvCenter->isMeshVertex()) {
        Vertex *v = plvCenter->getMeshVertex();
        foreach (int eId, *v->getIncidentEdgeIds()) {
            Edge *e = m_originalMesh->getEdge(eId);
            if (meshVertices.contains(e->getConnectedVertex(v)) && !meshEdges.contains(e))
                meshEdges.append(e);
        }
    } else if (plvCenter->isCrossingVertex()) {
        Edge *e = plvCenter->getCrossingEdge();
        if (meshVertices.contains(e->getVertex(0)) && meshVertices.contains(e->getVertex(1)) && !meshEdges.contains(e))
            meshEdges.append(e);
    }
    foreach (CellEdge *ce, incidentEdges) {
        const int numberOfPolyVertices = ce->getPolylineVertices()->size();
        for (int i = 1; i < numberOfPolyVertices - 1; i++) {
            PolylineVertex *plv = ce->getPolylineVertices()->at(i);
            if (plv->isMeshVertex()) {
                Vertex *v = plv->getMeshVertex();
                foreach (int eId, *v->getIncidentEdgeIds()) {
                    Edge *e = m_originalMesh->getEdge(eId);
                    if (meshVertices.contains(e->getConnectedVertex(v)) && !meshEdges.contains(e))
                        meshEdges.append(e);
                }
            } else if (plv->isCrossingVertex()) {
                Edge *e = plv->getCrossingEdge();
                if (meshVertices.contains(e->getVertex(0)) && meshVertices.contains(e->getVertex(1)) &&!meshEdges.contains(e))
                    meshEdges.append(e);
            }
        }
    }

    //Create new cell face
    int newFaceId = incidentFaceIds->at(0);
    CellFace *cfNew = new CellFace(borderCellEdges, meshEdges, meshVertices, newFaceId);
    m_cellFaces[newFaceId] = cfNew;

    //Update incident face ids of border edges and vertices
    const int numberOfBorderEdges = borderCellEdges.size();
    for (int i = 0; i < numberOfBorderEdges; i++) {
        const int oldId = faceIdsForUpdate[i];
        const int nextOldId = faceIdsForUpdate[(i + 1) % numberOfBorderEdges];
        ((CellEdge *) borderCellEdges.at(i))->updateIncidentFaceId(oldId, incidentFaceIds->at(0));
        CellVertex *cvShared = (CellVertex *) borderCellEdges[i]->getSharedVertex(borderCellEdges[(i+1)%numberOfBorderEdges]);
        if (oldId == nextOldId) {
            cvShared->updateIncidentFaceId(oldId, newFaceId);
        } else {
            cvShared->updateIncidentFaceId(nextOldId, -1);
            cvShared->updateIncidentFaceId(oldId, newFaceId);
        }
    }

    //Update incident edge ids of border vertices
    for (int i = 0; i < numberOfIncidentEdges; i++) {
        CellVertex *cv = (CellVertex *) incidentEdges.at(i)->getConnectedVertex(centerVertex);
        cv->updateIncidentEdgeId(incidentEdges.at(i)->getId(), -1);
    }


    //Delete old faces and free up their indices
    for (int i = 1; i < numberOfIncidentFaces; i++) {   //Skip the first ID
        const int faceId = incidentFaceIds->at(i);
        m_freeFaceIds.append(faceId);
        if (m_cellFaceParameterizations[faceId])
            delete m_cellFaceParameterizations[faceId];
        m_cellFaceParameterizations[faceId] = 0;
        m_cellFaces[faceId] = 0;
        delete incidentFaces[i];
    }

    if (m_cellFaceParameterizations[newFaceId])
        delete m_cellFaceParameterizations[newFaceId];
    m_cellFaceParameterizations[newFaceId] = 0;
    delete incidentFaces[0];

    //Delete old edges and free up their indices
    for (int i = 0; i < numberOfIncidentEdges; i++) {
        const int edgeId = incidentEdges.at(i)->getId();
        m_freeEdgeIds.append(edgeId);
        m_cellEdges[edgeId] = 0;
        delete incidentEdges[i];
    }

    m_freeVertexIds.append(centerVertex->getId());
    m_cellVertices[centerVertex->getId()] = 0;
    delete centerVertex;

    return cfNew;
}

CellEdge *CellMesh::mergeEdgesAroundUnecessaryVertex(CellVertex *centerVertex)
{
    if (centerVertex->getIncidentEdgeIds()->size() != 2 || centerVertex->getIncidentFaceIds()->size() != 2)
        return 0;

    const int cvId = centerVertex->getId();

    const int ce0Id = centerVertex->getIncidentEdgeIds()->at(0);    //will be used as new ID
    const int ce1Id = centerVertex->getIncidentEdgeIds()->at(1);    //will be deleted

    const int cf0Id = centerVertex->getIncidentFaceIds()->at(0);
    const int cf1Id = centerVertex->getIncidentFaceIds()->at(1);

    CellEdge *ce0 = m_cellEdges[ce0Id];
    CellEdge *ce1 = m_cellEdges[ce1Id];

    CellVertex *cv0 = (CellVertex *) ce0->getConnectedVertex(centerVertex);
    CellVertex *cv1 = (CellVertex *) ce1->getConnectedVertex(centerVertex);

    CellFace *cf0 = m_cellFaces[cf0Id];
    CellFace *cf1 = m_cellFaces[cf1Id];

    const int numPLV0 = ce0->getPolylineVertices()->size();
    const int numPLV1 = ce1->getPolylineVertices()->size();
    QVector<PolylineVertex *> PLVs;
    PLVs.reserve(numPLV0 + numPLV1 - 1);  //only allocates memory -> no resize
    if ((CellVertex *) ce0->getVertex(1) == centerVertex)
        for (int i = 0; i < numPLV0 - 1; i++)
            PLVs.append(new PolylineVertex(ce0->getPolylineVertices()->at(i)));
    else
        for (int i = ce0->getPolylineVertices()->size() - 1; i >= 1 ; i--)
            PLVs.append(new PolylineVertex(ce0->getPolylineVertices()->at(i)));
    //note that we do not add the last PLV since a copy of it is also contained in the next celledge

    if ((CellVertex *) ce1->getVertex(0) == centerVertex)
        for (int i = 0; i < numPLV1; i++)
            PLVs.append(new PolylineVertex(ce1->getPolylineVertices()->at(i)));
    else
        for (int i = ce1->getPolylineVertices()->size() - 1; i >= 0 ; i--)
            PLVs.append(new PolylineVertex(ce1->getPolylineVertices()->at(i)));

    CellEdge *ceNew = new CellEdge(cv0, cv1, PLVs, ce0Id);
    m_cellEdges[ce0Id] = ceNew;

    ceNew->addIncidentFaceId(cf0Id);
    ceNew->addIncidentFaceId(cf1Id);

    cf0->mergeEdges(ce0, ce1, ceNew);
    cf1->mergeEdges(ce0, ce1, ceNew);

    cv1->updateIncidentEdgeId(ce1Id, ce0Id);

    m_cellFaceParameterizations[cf0Id] = 0;
    m_cellFaceParameterizations[cf1Id] = 0;

    m_freeEdgeIds.append(ce1Id);
    m_cellEdges[ce1Id] = 0;
    delete ce1;

    m_freeVertexIds.append(cvId);
    m_cellVertices[cvId] = 0;
    delete centerVertex;

    return ceNew;
}

bool CellMesh::checkForConsistentOrientation()
{
    bool allOk = true;
    foreach (CellEdge *ce, m_cellEdges) {
        CellFace *cf0 = m_cellFaces[ce->getIncidentFaceId(0)];
        CellFace *cf1 = m_cellFaces[ce->getIncidentFaceId(1)];
        allOk = allOk && (cf0->edgeIsInverted(cf0->getInternalEdgeId(ce)) != cf1->edgeIsInverted(cf1->getInternalEdgeId(ce)));
    }

    return allOk;
}

void CellMesh::checkAndEnforceCellOrientation()
{
    bool allOk = this->checkForConsistentOrientation();

    if (!allOk) {
        foreach (CellFace *cf, m_cellFaces) {
            const int numberOfCellEdges = cf->getNumberOfEdges();
            //This method is not entierely stable, so check for each cell vertex and count right/wrong orientations
            int correctOrientation = 0;
            int wrongOrientation = 0;
            for (int iCe = 0; iCe < numberOfCellEdges; iCe++) {
                CellEdge *ce0 = (CellEdge *) cf->getEdge(iCe);
                CellEdge *ce1 = (CellEdge *) cf->getEdge((iCe + 1) % numberOfCellEdges);
                PolylineVertex *plvCorner = ce1->getPolylineVertices()->first();
                PolylineVertex *plv0 = ce0->getPolylineVertices()->at(ce0->getPolylineVertices()->size() - 2);
                PolylineVertex *plv1 = ce1->getPolylineVertices()->at(1);

                if (cf->edgeIsInverted(0))
                    plv0 = ce0->getPolylineVertices()->at(1);
                if (cf->edgeIsInverted(1)) {
                    plv1 = ce1->getPolylineVertices()->at(ce1->getPolylineVertices()->size() - 2);
                    plvCorner = ce1->getPolylineVertices()->last();
                }

                QVector3D cornerNormal;
                const int numIncidentFaces = plvCorner->getIncidentFaceIds()->size();
                for (int i = 0; i < numIncidentFaces; i++)
                    cornerNormal += m_originalMesh->getFace(plvCorner->getIncidentFaceIds()->at(i))->getNormal();

                QVector3D edgeNormal = QVector3D::crossProduct(plvCorner->getPosition() - plv0->getPosition(), plv1->getPosition() - plvCorner->getPosition());

                if (QVector3D::dotProduct(edgeNormal, cornerNormal) < 0)
                    wrongOrientation++;
                else
                    correctOrientation++;

            }
            if (wrongOrientation > correctOrientation) {
                cf->switchOrientation();
                m_cellFaceParameterizations[cf->getId()] = 0;
            }
        }
        enforceConsistentOrientation(); //Since the method is still not stable, make at least sure that the orientation is consistent
    }

}

void CellMesh::enforceConsistentOrientation()
{
    const int numberOfCellFaces = m_cellFaces.size();
    QVector<bool> checked(numberOfCellFaces, false);

    QList<CellFace *> queue;
    queue.append(m_cellFaces[0]);
    checked[0] = true;

    while (!queue.isEmpty()) {
        CellFace *cf = queue.takeFirst();
        const int numCfEdges = cf->getNumberOfEdges();
        for (int i = 0; i < numCfEdges; i++) {
            CellEdge *ce = (CellEdge *) cf->getEdge(i);
            CellFace *cfNeighbor = m_cellFaces[cf->getEdge(i)->getConnectedFaceId(cf->getId())];
            const int cfNeighborId = cfNeighbor->getId();
            if (!checked[cfNeighborId]) {
                queue.append(cfNeighbor);
                checked[cfNeighborId] = true;
            }
            if (cf->edgeIsInverted(i) == cfNeighbor->edgeIsInverted(cfNeighbor->getInternalEdgeId(ce))) {
                cfNeighbor->switchOrientation();
                m_cellFaceParameterizations[cfNeighborId] = 0;
            }
        }

        if (queue.isEmpty()) {
            int newConnectedComp = checked.indexOf(false);
            if (newConnectedComp != -1) {
                queue.append(m_cellFaces[newConnectedComp]);
                checked[newConnectedComp] = true;
            }
        }
    }

}

void CellMesh::invertTotalOrientation()
{
    this->cleanUpMeshIndices();

    const int numberOfCellFaces = m_cellFaces.size();
    for (int i = 0; i < numberOfCellFaces; i++) {
        m_cellFaces[i]->switchOrientation();
        m_cellFaceParameterizations[i] = 0;
    }

    this->addInfo(QString("Orientation of cellfaces inverted"));
}

QVector<int> CellMesh::getValencyHistogramCellVertices()
{
    QVector<int> histogram;
    foreach(CellVertex *cv, m_cellVertices) {
        if (cv) {
            const int val = cv->getIncidentEdgeIds()->size();
            if (histogram.size() < val + 1)
                histogram.resize(val + 1);
            histogram[val] += 1;
        }
    }
    return histogram;
}

QVector<int> CellMesh::getValencyHistogramCellFaces()
{
    QVector<int> histogram;
    foreach(CellFace *cf, m_cellFaces) {
        if (cf) {
            const int val = cf->getNumberOfEdges();
            if (histogram.size() < val + 1)
                histogram.resize(val + 1);
            histogram[val] += 1;
        }
    }
    return histogram;
}

//void CellMesh::subdivideSpheres()
//{
//    //TODO this procedure currently takes only one sphere -> detect all spheres
//    CellFace *firstHalf = 0;
//    CellFace *secondHalf = 0;
//    for (int i = 0; i < m_cellFaces.size(); i++) {
//        if (m_cellFaces[i]->getNumberOfVertices() == 1) {
//            firstHalf = m_cellFaces[i];
//            secondHalf = this->getNeighbor(firstHalf, 0);
//            break;
//        }
//    }

//    if (firstHalf == 0 || secondHalf == 0)
//        return;

//    //TODO parameterize both cells

//    //TODO add new corner points (four inside and four on the border)

//    //TODO build subdivided cells

//    //TODO update cell structure

//    //TODO DEBUG test if subdivided cells are still "correct" (i.e. connected components)
//}

int CellMesh::getNumberOfVertices()
{
    return m_cellVertices.size() - m_cellVertices.count(0);
}

CellVertex *CellMesh::getVertex(int id)
{
    return m_cellVertices[id];
}

int CellMesh::getNumberOfEdges()
{
    return m_cellEdges.size() - m_cellEdges.count(0);
}

CellEdge *CellMesh::getEdge(int id)
{
    return m_cellEdges[id];
}

int CellMesh::getNumberOfFaces()
{
    return m_cellFaces.size() - m_cellFaces.count(0);
}

CellFace *CellMesh::getFace(int id)
{
    return m_cellFaces[id];
}

int CellMesh::getNeighborId(CellFace *face, int number)
{
    Edge *e = face->getEdge(number);
    int id = e->getIncidentFaceId(0);
    if (id == face->getId())
        id = e->getIncidentFaceId(1);

    return id;
}

CellFace *CellMesh::getNeighbor(CellFace *face, int number)
{
    Edge *e = face->getEdge(number);
    int id = e->getIncidentFaceId(0);
    if (id == face->getId())
        id = e->getIncidentFaceId(1);

    if (id == -1)
        return 0;
    else
        return m_cellFaces[id];
}

void CellMesh::setParameterInterval(QPair<double, double> parameterInterval)
{
    m_parameterInterval = parameterInterval;
}

QPair<double, double> CellMesh::getParameterInterval()
{
    return m_parameterInterval;
}

Parameterization *CellMesh::getCellFaceParameterization(int i)
{
    if (!m_cellFaces[i])
        return 0;

    if (!m_cellFaceParameterizations[i])
        m_cellFaceParameterizations[i] = Parameterization::harmonicMap(m_cellFaces[i], m_originalMesh, m_parameterizationWeightType, m_parameterizationBorderType, m_parameterInterval, m_parameterizationStorageReservation_ExpectedEdgesPerVertex);

    return m_cellFaceParameterizations[i];
}

Parameterization *CellMesh::getQuadCellParameterization(int i)
{
    if (!m_cellFaces[i])
        return 0;

    if (m_cellFaces[i] && (!m_cellFaceParameterizations[i] ||
                           (m_cellFaceParameterizations[i]->getBorderType() != Parameterization::BT_equidistantUnitSquare &&
                            m_cellFaceParameterizations[i]->getBorderType() != Parameterization::BT_lengthWeightedUnitSquare)))
        m_cellFaceParameterizations[i] = Parameterization::harmonicMap(m_cellFaces[i], m_originalMesh, m_parameterizationWeightType, Parameterization::BT_lengthWeightedUnitSquare, m_parameterInterval, m_parameterizationStorageReservation_ExpectedEdgesPerVertex);

    return m_cellFaceParameterizations[i];
}

void CellMesh::updateBoundaryParameterization(int ceId, const QVector<double> &newParam, const bool onlyInteriorPoints)
{
    CellEdge * const ce = m_cellEdges[ceId];
    const int numPLV = ce->getPolylineVertices()->size();
    if (numPLV != newParam.size())
        return;

    for (int k = 0; k < 2; k++) {
        const int incidentFaceId = ce->getIncidentFaceId(k);
        Parameterization * const incidentCfParam = this->getQuadCellParameterization(incidentFaceId);
        const QVector2D startParam = incidentCfParam->getParameter(ce->getPolylineVertices()->first());
        const QVector2D endParam = incidentCfParam->getParameter(ce->getPolylineVertices()->last());
        const QVector2D borderVectorNormalized = (endParam - startParam).normalized();

        int skipPLV = 0;
        if (onlyInteriorPoints)
            skipPLV = 1;

        for (int i = skipPLV; i < numPLV - skipPLV; i++) {
            const QVector2D param2D = startParam + newParam[i] * borderVectorNormalized;
            incidentCfParam->updateParameter(ce->getPolylineVertices()->at(i), param2D.x(), param2D.y());
        }
    }
}

void CellMesh::setParameterizationWeightType(Parameterization::WeightType wt)
{
    m_parameterizationWeightType = wt;
}

void CellMesh::setParameterizationBorderType(Parameterization::BorderType bt)
{
    m_parameterizationBorderType = bt;
}

int CellMesh::getParameterizationStorageReservation()
{
    return m_parameterizationStorageReservation_ExpectedEdgesPerVertex;
}

void CellMesh::setParameterizationStorageReservation(int expectedEdgesPerVertex)
{
    if (expectedEdgesPerVertex > 0)
        m_parameterizationStorageReservation_ExpectedEdgesPerVertex = expectedEdgesPerVertex;
}

double CellMesh::limitQuadCellParameterizationsToParameterIntervall()
{
    this->cleanUpMeshIndices();
    const int numberOfCellFaces = m_cellFaces.size();
    double totalDif = 0;
    for (int i = 0; i < numberOfCellFaces; i++)
        totalDif += this->getQuadCellParameterization(i)->limitToParameterInterval();

    return totalDif;
}

QVector<double> CellMesh::calculateApproximateSurfaceAreas()
{
    qDebug() << "Calculation of surface areas for individual cells is not accurate, if cell edges contain mesh vertices!";
    this->cleanUpMeshIndices();

    const int numberOfCellFaces = m_cellFaces.size();
    QVector<double> areas(numberOfCellFaces);

    for (int cfId = 0; cfId < numberOfCellFaces; cfId++) {
        CellFace *cf = m_cellFaces[cfId];
        const int numMeshVertices = cf->getMeshVertices()->size();
        double cfArea = 0;

        //check mesh vertices
        for (int k = 0; k < numMeshVertices; k++) {
            Vertex *v = cf->getMeshVertices()->at(k);
            const int numIncidentFaces = v->getIncidentFaceIds()->size();
            for (int i = 0; i < numIncidentFaces; i++) {
                Face *f = m_originalMesh->getFace(v->getIncidentFaceIds()->at(i));
                cfArea += f->calculateAreaByCorners() / f->getNumberOfVertices();   //that way the area is counted fully, if all incident vertices are in the same cell
            }
        }


        //check potential mesh vertices along cell edge
        const int numIncidentCellEdges = cf->getNumberOfEdges();
        for (int iCe = 0; iCe < numIncidentCellEdges; iCe++) {
            CellEdge *ce = (CellEdge *) cf->getEdge(iCe);

            const int numPLVs = ce->getPolylineVertices()->size();

            for (int k = 0; k < numPLVs; k++) {
                PolylineVertex *plv = ce->getPolylineVertices()->at(k);

                if (plv->isMeshVertex()) {  //Only consider mesh vertices
                    Vertex *v = plv->getMeshVertex();
                    const int numIncidentFaces = v->getIncidentFaceIds()->size();
                    for (int i = 0; i < numIncidentFaces; i++) {
                        Face *f = m_originalMesh->getFace(v->getIncidentFaceIds()->at(i));
                        const int numFaceVertices = f->getNumberOfVertices();
                        for (int j = 0; j < numFaceVertices; j++) {
                            Vertex *fv = f->getVertex(j);
                            //TODO that is not 100% accurate
                            //...will produce wrong results if two mesh verices are on the same mesh edge, but no vertex is inside the cell
                            if (fv != v && this->getCellFaceParameterization(cfId)->contains(fv)) { //Only add area if at least one vertex of the face is in the cell
                                cfArea += f->calculateAreaByCorners() / numFaceVertices;
                                break;
                            }
                        }
                    }
                }

            }
        }


        areas[cfId] = cfArea;
    }

    return areas;
}

void CellMesh::translate(const QVector3D &translationVector)
{
    this->cleanUpMeshIndices();

    foreach (CellVertex *cv, m_cellVertices)
        cv->updatePosition(cv->getPosition() + translationVector);

    foreach (CellEdge *ce, m_cellEdges)
        foreach (PolylineVertex *plv, *ce->getPolylineVertices())
            plv->updatePosition(plv->getPosition() + translationVector);

}

void CellMesh::scale(double factorX, double factorY, double factorZ)
{
    this->cleanUpMeshIndices();

    foreach (CellVertex *cv, m_cellVertices) {
        const QVector3D pos = cv->getPosition();
        cv->updatePosition(QVector3D(pos.x() * factorX, pos.y() * factorY, pos.z() * factorZ));
    }

    foreach (CellEdge *ce, m_cellEdges)
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            const QVector3D pos = plv->getPosition();
            plv->updatePosition(QVector3D(pos.x() * factorX, pos.y() * factorY, pos.z() * factorZ));
        }
}

Mesh *CellMesh::getOriginalMesh()
{
    return m_originalMesh;
}

QString CellMesh::getInfo()
{
    return m_info;
}

void CellMesh::resetInfo()
{
    m_info = QString();
}

void CellMesh::addInfo(QString additionalInfo)
{
    m_info.append(additionalInfo);
}

void CellMesh::cleanUpMeshIndices()
{
    //Clean up vertices
    m_freeVertexIds.clear();
    while (!m_cellVertices.last())
        m_cellVertices.remove(m_cellVertices.size() - 1);

    int index = m_cellVertices.indexOf(0);
    while (index != -1) {
        CellVertex *cv = m_cellVertices.last();
        cv->updateId(index);
        m_cellVertices[index] = cv;
        m_cellVertices.remove(m_cellVertices.size() - 1);
        while (!m_cellVertices.last())
            m_cellVertices.remove(m_cellVertices.size() - 1);

        index = m_cellVertices.indexOf(0);
    }

    //Clean up edges
    m_freeEdgeIds.clear();
    while (!m_cellEdges.last())
        m_cellEdges.remove(m_cellEdges.size() - 1);

    index = m_cellEdges.indexOf(0);
    while (index != -1) {
        CellEdge *ce = m_cellEdges.last();

        for (int i = 0; i < 2; i++)
            ((CellVertex *) ce->getVertex(i))->updateIncidentEdgeId(ce->getId(), index);

        ce->updateId(index);
        m_cellEdges[index] = ce;
        m_cellEdges.remove(m_cellEdges.size() - 1);
        while (!m_cellEdges.last())
            m_cellEdges.remove(m_cellEdges.size() - 1);

        index = m_cellEdges.indexOf(0);
    }

    //Clean up faces
    m_freeFaceIds.clear();
    while (!m_cellFaces.last()) {
        m_cellFaces.remove(m_cellFaces.size() - 1);
        m_cellFaceParameterizations.remove(m_cellFaces.size() - 1);
    }

    index = m_cellFaces.indexOf(0);
    while (index != -1) {
        CellFace *cf = m_cellFaces.last();

        const int numEdges = cf->getNumberOfEdges();
        for (int i = 0; i < numEdges; i++)
            ((CellEdge *) cf->getEdge(i))->updateIncidentFaceId(cf->getId(), index);

        const int numVertices = cf->getNumberOfVertices();
        for (int i = 0; i < numVertices; i++)
            ((CellVertex *) cf->getVertex(i))->updateIncidentFaceId(cf->getId(), index);

        cf->updateId(index);
        m_cellFaces[index] = cf;
        m_cellFaceParameterizations[index] = m_cellFaceParameterizations.last();
        m_cellFaces.remove(m_cellFaces.size() - 1);
        m_cellFaceParameterizations.remove(m_cellFaceParameterizations.size() - 1);
        while (!m_cellFaces.last())
            m_cellFaces.remove(m_cellFaces.size() - 1);

        index = m_cellFaces.indexOf(0);
    }
}

void CellMesh::checkMesh()
{
    std::cout << "MESH CHECK: "
              << "Number of faces: " << this->getNumberOfFaces()
              << " / Number of edges: " << this->getNumberOfEdges()
              << " / Number of vertices: " << this->getNumberOfVertices()
              << std::endl;

    for(int i = 0; i < m_cellFaces.size(); i++) {
        CellFace *cf = m_cellFaces[i];
        if (cf) {
            const int cfId = cf->getId();
            const int numVertices = cf->getNumberOfVertices();
            const int numEdges = cf->getNumberOfEdges();
            if (numVertices != numEdges)
                std::cout << "Face " << cfId << " has " << numVertices << " vertices and " << numEdges << " edges!" << std::endl;
            for (int j = 0; j < numVertices; j++)
                if (!cf->getVertex(j)->getIncidentFaceIds()->contains(cfId))
                    std::cout << "Face " << cfId << " knows vertex " << cf->getVertex(j)->getId() << ", but vertex does not know face" << std::endl;
            for (int j = 0; j < numEdges; j++)
                if (cf->getEdge(j)->getIncidentFaceId(0) != cfId && cf->getEdge(j)->getIncidentFaceId(1) != cfId)
                    std::cout << "Face " << cfId << " knows edge " << cf->getEdge(j)->getId() << ", but edge does only know face" << std::endl;
        }
    }

    for (int i = 0; i < m_cellEdges.size(); i++) {
        CellEdge *ce = m_cellEdges[i];
        if (ce) {
            const int ceId = ce->getId();
            for (int j = 0; j < 2; j++) {
                if (!ce->getVertex(j)->getIncidentEdgeIds()->contains(ceId))
                    std::cout << "Edge " << ceId << " knows vertex " << ce->getVertex(j)->getId() << ", but vertex does not know edge" << std::endl;
                CellFace *incidentFace = m_cellFaces[ce->getIncidentFaceId(j)];
                if (!incidentFace) {
                    std::cout << "Edge " << ceId << " is incident to face " << ce->getIncidentFaceId(j) << ", but face does not exist!" << std::endl;
                } else {
                    bool faceContainsEdge = false;
                    const int numEdges = incidentFace->getNumberOfEdges();
                    for (int k = 0; k < numEdges; k++)
                        faceContainsEdge = (faceContainsEdge || incidentFace->getEdge(k) == ce);
                    if (!faceContainsEdge)
                        std::cout << "Edge " << ceId << " knows face " << incidentFace->getId() << ", but faces does not contain edge" << std::endl;
                }
            }
        }
    }

    for (int i = 0; i < m_cellVertices.size(); i++) {
        CellVertex *cv = m_cellVertices[i];
        if (cv) {
            const int cvId = cv->getId();
            const int numEdges = cv->getIncidentEdgeIds()->size();
            const int numFaces = cv->getIncidentFaceIds()->size();

            for (int j = 0; j < numEdges; j++) {
                CellEdge *ce = m_cellEdges[cv->getIncidentEdgeIds()->at(j)];
                if (!ce)
                    std::cout << "Vertex " << cvId << " is incident to edge " << cv->getIncidentEdgeIds()->at(j) << ", but edge does not exist!" << std::endl;
                else if (ce->getVertex(0) != cv && ce->getVertex(1) != cv)
                    std::cout << "Vertex " << cvId << " knows edge " << ce->getId() << ", but edge does not contain vertex!" << std::endl;
            }

            for (int j = 0; j < numFaces; j++) {
                CellFace *cf = m_cellFaces[cv->getIncidentFaceIds()->at(j)];
                if (!cf) {
                    std::cout << "Vertex " << cvId << " is incident to face " << cv->getIncidentFaceIds()->at(j) << ", but face does not exist!" << std::endl;
                } else {
                    bool faceContainsVertex = false;
                    const int numVertices = cf->getNumberOfVertices();
                    for (int k = 0; k < numVertices; k++)
                        faceContainsVertex = (faceContainsVertex || cf->getVertex(k) == cv);
                    if (!faceContainsVertex)
                        std::cout << "Vertex " << cvId << " knows face " << cf->getId() << ", but face does not contain vertex!" << std::endl;
                }
            }
        }
    }
}

void CellMesh::saveToFile(QString filename, QString optionalHeader)
{
    this->cleanUpMeshIndices(); //needs to be done to avoid problems when loading the cell mesh

    const int numCellVertices = m_cellVertices.size();
    const int numCellEdges = m_cellEdges.size();
    const int numCellFaces = m_cellFaces.size();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    out << m_info;
    if (!m_info.endsWith("\n"))
        out << "\n";
    out << optionalHeader;
    if (!optionalHeader.endsWith("\n"))
        out << "\n";
    out << "===\n";
    out << QString::number(numCellVertices) << " " << QString::number(numCellEdges) << " " << QString::number(numCellFaces) << " "
        << QString::number(m_parameterizationWeightType) << " " << QString::number(m_parameterizationBorderType) << "\n"
        << QString::number(m_parameterInterval.first) << " " << QString::number(m_parameterInterval.second) << "\n";
    out << "===\n";
    for (int i = 0; i < numCellVertices; i++) {
        CellVertex *cv = m_cellVertices[i];
        out << QString::number(cv->getId()) << " " <<
               QString::number(cv->getPosX()) << " " <<
               QString::number(cv->getPosY()) << " " <<
               QString::number(cv->getPosZ()) << "\n";
    }

    out << "===\n";

    for (int i = 0; i < numCellEdges; i++) {
        CellEdge *ce = m_cellEdges[i];
        const int numPolyVert = ce->getPolylineVertices()->size();
        int isG1 = 1;
        if (!ce->isG1Edge())
            isG1 = 0;
        out << QString::number(ce->getId()) << " " <<
               QString::number(ce->getVertex(0)->getId()) << " " <<
               QString::number(ce->getVertex(1)->getId()) << " " <<
               QString::number(numPolyVert) << " " <<
               QString::number(isG1) << "\n";
        for (int j = 0; j < numPolyVert; j++) {
            PolylineVertex *plv = ce->getPolylineVertices()->at(j);
            int vId = -1;
            int eId = -1;
            int fId = -1;
            if (plv->isMeshVertex())
                vId = plv->getMeshVertex()->getId();
            else if (plv->isCrossingVertex())
                eId = plv->getCrossingEdge()->getId();
            else if (plv->isFreeVertex()){
                fId = plv->getIncidentFaceIds()->at(0);
            } else {
                qDebug() << "Error in cellmesh saving: Polyline vertex is neither mesh vertex, crossing vertex or free vertex!";
            }
            out << QString::number(plv->getPosX()) << " " <<
                   QString::number(plv->getPosY()) << " " <<
                   QString::number(plv->getPosZ()) << " " <<
                   QString::number(vId) << " " <<
                   QString::number(eId) << " " <<
                   QString::number(fId) << "\n";
        }
    }

    out << "===\n";

    for (int i = 0; i < numCellFaces; i++) {
        CellFace *cf = m_cellFaces[i];
        const int numCellFaceVert = cf->getNumberOfVertices();
        const int numContainedVertices = cf->getMeshVertices()->size();
        const int numContainedEdges = cf->getMeshEdges()->size();

        out << QString::number(cf->getId()) << " "
            << QString::number(numCellFaceVert);
        for (int j = 0; j < numCellFaceVert; j++)
            out << " " << QString::number(cf->getVertex(j)->getId());

        for (int j = 0; j < numCellFaceVert; j++)
            out << " " << QString::number(cf->getEdge(j)->getId());

        out << "\n" << QString::number(numContainedVertices);
        for (int j = 0; j < numContainedVertices; j++)
            out << " " << QString::number(cf->getMeshVertices()->at(j)->getId());

        out << "\n" << QString::number(numContainedEdges);
        for (int j = 0; j < numContainedEdges; j++)
            out << " " << QString::number(cf->getMeshEdges()->at(j)->getId());

        out << "\n";
    }

    out << "===";

    file.close();
}

CellMesh *CellMesh::loadFromFile(QString cellmeshFile, Mesh *origMesh)
{
    QFile file(cellmeshFile);
    if (!origMesh || !file.open(QIODevice::ReadOnly))
        return 0;
    QTextStream in(&file);

    CellMesh *cm = new CellMesh();

    QString info;
    QString line = in.readLine();
    //skip until first seperator
    while (line != "===") {
        info.append(line).append("\n");
        line = in.readLine();
    }
    cm->addInfo(info);

    //Read mesh dimensions
    line = in.readLine();   //This line contains the size information of the mesh
    QStringList split = line.split(" ");
    const int numCellVertices = split[0].toInt();
    const int numCellEdges = split[1].toInt();
    const int numCellFaces = split[2].toInt();
    cm->m_cellVertices = QVector<CellVertex *>(numCellVertices, 0);
    cm->m_cellEdges = QVector<CellEdge *>(numCellEdges, 0);
    cm->m_cellFaces = QVector<CellFace *>(numCellFaces, 0);
    cm->m_parameterizationWeightType = (Parameterization::WeightType) split[3].toInt();
    cm->m_parameterizationBorderType = (Parameterization::BorderType) split[4].toInt();
    cm->m_cellFaceParameterizations = QVector<Parameterization *>(numCellFaces, 0);
    cm->m_originalMesh = origMesh;
    cm->m_parameterizationStorageReservation_ExpectedEdgesPerVertex = origMesh->getMaxNumberOfEdgesOnVertex() + 1;

    line = in.readLine();   //next separator (old version) or parameter interval
    if (line != "===") {
        split = line.split(" ");
        const double paramA = split[0].toDouble();
        const double paramB = split[1].toDouble();
        cm->m_parameterInterval = QPair<double, double>(paramA, paramB);
        line = in.readLine();   //next separator
    } else {
        cm->m_parameterInterval = QPair<double, double>(0, 1);   //default value for old files
    }

    //read cell vertices
    line = in.readLine();   //read the first vertex
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        QVector3D pos = QVector3D(split[1].toFloat(), split[2].toFloat(), split[3].toFloat());
        CellVertex *cv = new CellVertex(pos, id);
        cm->m_cellVertices[id] = cv;

        line = in.readLine();
    }

    //read cell edges
    line = in.readLine();   //read the first edge
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        const int cv0Id = split[1].toInt();
        const int cv1Id = split[2].toInt();
        CellVertex *cv0 = cm->m_cellVertices[cv0Id];
        CellVertex *cv1 = cm->m_cellVertices[cv1Id];
        cv0->addIncidentEdgeId(id);
        cv1->addIncidentEdgeId(id);

        const int numPLV = split[3].toInt();
        QVector<PolylineVertex *> polyVertices(numPLV);
        bool isG1 = true;
        if (split.size() > 4 && split[4].toInt() == 0)  //the first condition ensures compability to older versions
            isG1 = false;

        for (int i = 0; i < numPLV; i++) {
            line = in.readLine();
            split = line.split(" ");
            QVector3D pos(split[0].toFloat(), split[1].toFloat(), split[2].toFloat());
            int vId = split[3].toInt();
            int eId = split[4].toInt();
            int fId = split[5].toInt();
            PolylineVertex *plv;
            if (vId != -1)
                plv = new PolylineVertex(origMesh->getVertex(vId));
            else if (eId != -1)
                plv = new PolylineVertex(pos, origMesh->getEdge(eId));
            else
                plv = new PolylineVertex(pos, fId);

            polyVertices[i] = plv;
        }

        CellEdge *ce = new CellEdge(cv0, cv1, polyVertices, id);
        ce->setG1Edge(isG1);

        cm->m_cellEdges[id] = ce;

        line = in.readLine();
    }

    //read cell faces
    line = in.readLine();   //read the first face
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        const int numCellFaceVert = split[1].toInt();
        QVector<Edge *> faceCellEdges(numCellFaceVert);
        for (int i = 0; i < numCellFaceVert; i++) {
            Edge *ce = cm->getEdge(split[2 + numCellFaceVert + i].toInt());
            faceCellEdges[i] = ce;

            cm->m_cellVertices[split[2 + i].toInt()]->addIncidentFaceId(id);    //add incident face to vertex
            ce->addIncidentFaceId(id);                                          //add incident face to edge
        }
        line = in.readLine();   //mesh vertices
        split = line.split(" ");
        const int numMeshVert = split[0].toInt();
        QVector<Vertex *> meshVertices(numMeshVert);
        for (int i = 0; i < numMeshVert; i++)
            meshVertices[i] = origMesh->getVertex(split[i+1].toInt());

        line = in.readLine();   //mesh edges
        split = line.split(" ");
        const int numMeshEdges = split[0].toInt();
        QVector<Edge *> meshEdges(numMeshEdges);
        for (int i = 0; i < numMeshEdges; i++)
            meshEdges[i] = origMesh->getEdge(split[i+1].toInt());

        cm->m_cellFaces[id] = new CellFace(faceCellEdges, meshEdges, meshVertices, id);

        line = in.readLine();   //next cell face
    }

    file.close();

    return cm;
}

void CellMesh::saveParameterizations(QString genericFilename, QString exportFilename, QString optionalHeader)
{
    const int numberOfCells = m_cellFaces.size();

    for (int i = 0; i < numberOfCells; i++) {
        this->getCellFaceParameterization(i)->saveToFile(genericFilename.arg(i), optionalHeader);

        if (!exportFilename.isEmpty())
            this->getCellFaceParameterization(i)->saveAsImage(exportFilename.arg(i));
    }
}

void CellMesh::loadParameterizations(QString genericFilename)
{
    const int numberOfCellFaces = m_cellFaces.size();
    for (int i = 0; i < numberOfCellFaces; i++)
        if (QFile(genericFilename.arg(i)).exists())
            m_cellFaceParameterizations[i] = Parameterization::loadFromFile(genericFilename.arg(i), m_cellFaces[i], m_originalMesh);
}

CellMesh::CellMesh()
{
    //Empty contructor -> to be used by the static load method
}

int CellMesh::takeFreeVertexId()
{
    if (m_freeVertexIds.isEmpty()) {
        m_cellVertices.append(0);
        return m_cellVertices.size()-1;
    } else {
        return m_freeVertexIds.takeFirst();
    }
}

int CellMesh::takeFreeEdgeId()
{
    if (m_freeEdgeIds.isEmpty()) {
        m_cellEdges.append(0);
        return m_cellEdges.size()-1;
    } else {
        return m_freeEdgeIds.takeFirst();
    }
}

int CellMesh::takeFreeFaceId()
{
    if (m_freeFaceIds.isEmpty()) {
        m_cellFaces.append(0);
        m_cellFaceParameterizations.append(0);
        return m_cellFaces.size()-1;
    } else {
        return m_freeFaceIds.takeFirst();
    }
}

CellVertex *CellMesh::subdivideEdge(CellEdge *ce, PolylineVertex *splitVertex, QVector<PolylineVertex *> *updateVertices)
{
    //Get parameterizations of incident cell faces to update the pointerToParameter-map
    QVector<Parameterization *> incidentParameterizations;
    for (int i = 0; i < 2; i++) {
        Parameterization *param = m_cellFaceParameterizations[ce->getIncidentFaceId(i)];
        if (param)
            incidentParameterizations.append(param);
    }
    const int existingParams = incidentParameterizations.size();

    //Split the list of mesh edges belonging to the cell edge
    QVector<PolylineVertex *> polyVertices0;
    QVector<PolylineVertex *> polyVertices1;
    QVector<PolylineVertex *> *polyVerticesFull = ce->getPolylineVertices();
    const int numberOfPolyVertices = polyVerticesFull->size();

    if (splitVertex == polyVerticesFull->first() || splitVertex == polyVerticesFull->last())
        return 0;

    int i = 0;
    PolylineVertex *plv = 0;
    for (; i < numberOfPolyVertices && splitVertex != plv; i++) {
        plv = polyVerticesFull->at(i);
        PolylineVertex *newPlv = new PolylineVertex(plv);
        polyVertices0.append(newPlv);
        if (updateVertices) {
            const int id = updateVertices->indexOf(plv);
            if (id > -1)
                (*updateVertices)[id] = newPlv;
        }
        for (int j = 0; j < existingParams; j++)
            incidentParameterizations[j]->addCopy(plv, newPlv);
    }

    PolylineVertex *newPlv = new PolylineVertex(plv);
    polyVertices1.append(newPlv);   //Needs to be the endpoint of both new cell edges
    for (int j = 0; j < existingParams; j++)
        incidentParameterizations[j]->addCopy(plv, newPlv);

    for (;i < numberOfPolyVertices; i++) {
        plv = polyVerticesFull->at(i);
        PolylineVertex *newPlv = new PolylineVertex(plv);
        polyVertices1.append(newPlv);
        if (updateVertices) {
            const int id = updateVertices->indexOf(plv);
            if (id > -1)
                (*updateVertices)[id] = newPlv;
        }
        for (int j = 0; j < existingParams; j++)
            incidentParameterizations[j]->addCopy(plv, newPlv);
    }

    //Create new cell vertex and add it to the list of cell vertices
    const int newVertexId = this->takeFreeVertexId();
    CellVertex *cvNew = new CellVertex(splitVertex->getPosition(), newVertexId);
    m_cellVertices[newVertexId] = cvNew;
    for (int j = 0; j < existingParams; j++)
        incidentParameterizations[j]->addCopy(splitVertex, cvNew);

    //Create new cell edges (one will replace the old edge, the other one will be added at the end of the list)
    const int originalEdgeId = ce->getId();
    const int newEdgeId = this->takeFreeEdgeId();
    CellVertex *cv0 = (CellVertex *) ce->getVertex(0);
    CellVertex *cv1 = (CellVertex *) ce->getVertex(1);
    CellEdge *split0 = new CellEdge(cv0, cvNew, polyVertices0, originalEdgeId);
    CellEdge *split1 = new CellEdge(cvNew, cv1, polyVertices1, newEdgeId);
    m_cellEdges[originalEdgeId] = split0;
    m_cellEdges[newEdgeId] = split1;

    //Set incident face IDs for new edges
    const int face0Id = ce->getIncidentFaceId(0);
    const int face1Id = ce->getIncidentFaceId(1);
    split0->addIncidentFaceId(face0Id);
    split0->addIncidentFaceId(face1Id);
    split1->addIncidentFaceId(face0Id);
    split1->addIncidentFaceId(face1Id);

    //Set the incident face IDs for the newly created cell vertex
    cvNew->addIncidentFaceId(face0Id);
    cvNew->addIncidentFaceId(face1Id);

    //Set/Update incident edge IDs for cell vertices (since the first edge has the same index as the old edge, the first cell vertex does not need to be updated)
    cv1->updateIncidentEdgeId(originalEdgeId, newEdgeId);
    cvNew->addIncidentEdgeId(originalEdgeId);
    cvNew->addIncidentEdgeId(newEdgeId);

    //Update the cell edges and cell vertices of the two incident faces
    m_cellFaces[face0Id]->replaceSplitEdges(ce, split0, split1);
    m_cellFaces[face1Id]->replaceSplitEdges(ce, split0, split1);

    //Remove the old edge
    delete ce;

    //Return newly created cell vertex
    return cvNew;
}

QVector3D CellMesh::calculateIntersectionWithLine(Vertex *v0, Vertex *v1, Parameterization *param, const QVector2D &projectedStart, const QVector2D &normal, QVector2D *intersectionParameter)
{
    QVector2D param0 = param->getParameter(v0);
    QVector2D param1 = param->getParameter(v1);
    //((lambda * v0 + (1 - lambda) * v1) - root) * normal = 0
    double lambda = QVector2D::dotProduct(projectedStart - param1, normal) / QVector2D::dotProduct(param0 - param1, normal);
    //To avoid numerical problems, set intersection to at least 2% of an edge length away from the vertices
    if (lambda < 0.02)
        lambda = 0.02;
    if (lambda > 0.98)
        lambda = 0.98;
    QVector3D intersection = lambda * v0->getPosition() + (1 - lambda) * v1->getPosition();
    (*intersectionParameter) = lambda * param0 + (1 - lambda) * param1;
    return intersection;
}

QVector2D CellMesh::calculateIntersectionCoefficients(const QVector2D &u0, const QVector2D &u1, const QVector2D &v0, const QVector2D &v1)
{
    const QVector2D u0u1 = u0 - u1;
    const QVector2D v1v0 = v1 - v0;
    const QVector2D v1u1 = v1 - u1;

    const double det = (u0u1.x() * v1v0.y()) - (u0u1.y() * v1v0.x());

    if (det == 0) {
        return QVector2D(-1, -1);
    } else {
        const double lambdaU = ((v1v0.y() * v1u1.x()) - (v1v0.x() * v1u1.y()))/det;
        const double lambdaV = ( -(u0u1.y() * v1u1.x()) + (u0u1.x() * v1u1.y()))/det;

        return QVector2D(lambdaU, lambdaV);
    }
}

void CellMesh::restrictToInterval(QVector2D *lambda, const float minX, const float maxX, const float minY, const float maxY)
{
    if (lambda->x() < minX)
        lambda->setX(minX);
    if (lambda->x() > maxX)
        lambda->setX(maxX);
    if (lambda->y() < minY)
        lambda->setY(minY);
    if (lambda->y() > maxY)
        lambda->setY(maxY);
}

int CellMesh::findSlice(double angle, QVector<double> &separationAngles)
{
    int index = 0;
    const int numberOfSlices = separationAngles.size();
    if (angle < separationAngles[0]) {
        index = numberOfSlices - 1;
    } else {
        while (index < numberOfSlices && angle >= separationAngles[index])
            index++;
        index--;
    }

    return index;
}

QVector<double> CellMesh::getBoundaryParameterizationFromIncidentQuadParam(int ceId)
{
   CellEdge * const ce = m_cellEdges[ceId];
   const int incidentCfId = ce->getIncidentFaceId(0);
   CellFace * const incidentCf = m_cellFaces[incidentCfId];
   return this->getQuadCellParameterization(incidentCfId)->getParameterizationOfBoundaryCurve(incidentCf->getInternalEdgeId(ce));
}

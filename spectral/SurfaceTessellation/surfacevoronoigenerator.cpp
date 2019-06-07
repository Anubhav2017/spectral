#include "surfacevoronoigenerator.h"

#include "priorityqueue.h"

SurfaceVoronoiGenerator::SurfaceVoronoiGenerator(Mesh *m, const bool initializeWithFirstFace):
    m_mesh(m), m_numberOfMeshFaces(m->getNumberOfFaces())
{
    this->initializeEmptyDijkstraVoronoi();

    if (initializeWithFirstFace)
        this->addSource(m_mesh->getFace(0));
}

void SurfaceVoronoiGenerator::initializeEmptyDijkstraVoronoi()
{
    const double infinity = std::numeric_limits<double>::infinity();
    m_voronoiCells.clear();
    m_voronoiId = QVector<int>(m_numberOfMeshFaces, -1);
    m_dijkstraDist = QVector<double>(m_numberOfMeshFaces, infinity);
    m_dijkstraPred = QVector<Face *>(m_numberOfMeshFaces, 0);
    m_sources.clear();
    m_priorityQueue = PriorityQueue<Face *>();
    m_inQueue = QVector<bool>(m_numberOfMeshFaces, false);
}

CellMesh *SurfaceVoronoiGenerator::generateCellMesh(bool atLeast3PlvPerEdge)
{
    return new CellMesh(m_mesh, m_voronoiId, m_sources.size(), atLeast3PlvPerEdge);
}

void SurfaceVoronoiGenerator::computeDiskLikeVoronoiTessellation(const int minNumberOfCells, bool instantTopologyCheck)
{
    bool allOk = false;
    while (!allOk) {
        allOk = true;

        while (!m_priorityQueue.isEmpty() && allOk) {
            Face *f = m_priorityQueue.dequeue();
            const int fId = f->getId();
            m_inQueue[fId] = false;
            const int numFaceEdges = f->getNumberOfEdges();

            const int predId = m_dijkstraPred[fId]->getId();
            const int fVoroId = m_voronoiId[predId];
            m_voronoiId[fId] = fVoroId;

            if (instantTopologyCheck && predId != fId) {    //topology check only if f is not a source
                //change in euler number equals new faces - new edges + new vertices
                int eulerNumberChange = 1;  //+1 face

                for (int i = 0; i < numFaceEdges; i++) {
                    Face *neighbor = m_mesh->getNeighbor(f, i);
                    if (m_voronoiId[neighbor->getId()] != fVoroId) //edge was not in cell
                        eulerNumberChange--;

                    Vertex *v = f->getVertex(i);
                    const int numIncidentFaces = v->getIncidentFaceIds()->size();
                    bool vertexWasInCell = false;
                    for (int k = 0; k < numIncidentFaces && !vertexWasInCell; k++) {
                        const int incidentFaceId = v->getIncidentFaceIds()->at(k);
                        if (fId != incidentFaceId && m_voronoiId[incidentFaceId] == fVoroId)
                            vertexWasInCell = true;
                    }

                    if (!vertexWasInCell)
                        eulerNumberChange++;
                }

                if (eulerNumberChange != 0) {
                    this->addSource(f);
                    allOk = false;
                }
            }

            if (allOk) {    //this can only be false if instant topology check is used
                for (int i = 0; i < numFaceEdges; i++) {
                    Face *neighbor = m_mesh->getNeighbor(f, i);
                    const int neighborId = neighbor->getId();
                    const double distToNeighbor = (f->getCenter() - neighbor->getCenter()).length();  //distance metric
                    const double dist = m_dijkstraDist[neighborId];
                    const double newDist = m_dijkstraDist[fId] + distToNeighbor;
                    if (newDist < dist) {
                        if (!m_inQueue[neighborId]) {
                            m_priorityQueue.enqueue(neighbor, newDist);
                            m_inQueue[neighborId] = true;
                        } else
                            m_priorityQueue.refreshDistance(neighbor, newDist);
                        m_dijkstraDist[neighborId] = newDist;
                        m_dijkstraPred[neighborId] = f;
                    }
                }
            }
        }

        //Dijkstra done, setup list of faces for topology check and further processing (unless there is still a connected component missing)
        const int numberOfVoronoiCells = m_sources.size();
        m_voronoiCells = QVector<QVector<Face *> >(numberOfVoronoiCells, QVector<Face *>());
        for (int i = 0; i < m_numberOfMeshFaces && allOk; i++) {
            const int voroId = m_voronoiId[i];
            Face *f = m_mesh->getFace(i);
            if (voroId != -1) {
                m_voronoiCells[voroId].append(f);
            } else {
                allOk = false;
                this->addSource(f);
            }
        }

        //Topology check
        if (!instantTopologyCheck) {    //only necessary if we did not use the instant topology check earlier
            for (int i = 0; i < numberOfVoronoiCells && allOk; i++) {
                const int eulerNumber = calculateEulerNumberOfCell(i);
                if (eulerNumber != 1) { //i.e. if topology is not like a disk
                    const int highestDistListID = findHighestDijkstraDistIndexInCell(i);
                    Face *newSource = m_voronoiCells[i][highestDistListID];

                    this->addSource(newSource);
                    allOk = false;
                }
            }
        }

        //Check if minimum number of cells is reached
        if (numberOfVoronoiCells < minNumberOfCells && allOk) {
            this->addSourceFromLargestCellWithHigestDistance();
            allOk = false;
        }
    }
}



bool SurfaceVoronoiGenerator::isSource(Face *f)
{
    return m_sources.contains(f);
}

void SurfaceVoronoiGenerator::addSource(Face *source)
{
    const int sourceId = source->getId();
    m_voronoiId[sourceId] = m_sources.size();
    m_sources.append(source);
    m_dijkstraPred[sourceId] = source;
    m_dijkstraDist[sourceId] = 0;

    m_priorityQueue.enqueue(source, 0);
    m_inQueue[sourceId] = true;
}

void SurfaceVoronoiGenerator::addSourceWithLargestDistanceInCell(const int cellId)
{
    const int index = findHighestDijkstraDistIndexInCell(cellId);
    addSource(m_voronoiCells[cellId][index]);
}

void SurfaceVoronoiGenerator::addSourceFromLargestCellWithHigestDistance()
{
    const int cellId = this->findLargestCell();
    const int index = findHighestDijkstraDistIndexInCell(cellId);
    addSource(m_voronoiCells[cellId][index]);
}

int SurfaceVoronoiGenerator::calculateEulerNumberOfCell(const int cellId)
{
    const int numFacesInCell = m_voronoiCells[cellId].size();

    int eulerNumber = numFacesInCell;   // + number of faces
    foreach (Face *f, m_voronoiCells[cellId]) {
        const int fId = f->getId();
        const int numNeighbors = f->getNumberOfEdges();
        for (int j = 0; j < numNeighbors; j++) {
            // - number of edges
            const int neighborId = m_mesh->getNeighbor(f, j)->getId();
            if (fId > neighborId || m_voronoiId[neighborId] != cellId) //fId > neighborId makes sure that each edge is only counted once
                eulerNumber--;

            // + number of vertices
            Vertex *v = f->getVertex(j);
            const int numIncidentFaces = v->getIncidentFaceIds()->size();
            bool countV = true;
            for (int k = 0; k < numIncidentFaces && countV; k++) {
                const int incidentFaceId = v->getIncidentFaceIds()->at(k);
                if (incidentFaceId > fId && m_voronoiId[incidentFaceId] == cellId)   //make sure to count each vertex only once
                    countV = false;
            }

            if (countV)
                eulerNumber++;
        }
    }

    return eulerNumber;
}

int SurfaceVoronoiGenerator::findLargestCell()
{
    int largestSize = m_voronoiCells[0].size();
    int largestId = 0;
    const int numCells = m_voronoiCells.size();
    for (int i = 1; i < numCells; i++) {
        if (m_voronoiCells[i].size() > largestSize) {
            largestSize = m_voronoiCells[i].size();
            largestId = i;
        }
    }

    return largestId;
}

int SurfaceVoronoiGenerator::findHighestDijkstraDistIndexInCell(const int cellId)
{
    const int cellSize = m_voronoiCells[cellId].size();

    double highestDist = m_dijkstraDist[m_voronoiCells[cellId][0]->getId()];
    int highestDistListID = 0;
    for (int i = 1; i < cellSize; i++) {
        const double dist = m_dijkstraDist[m_voronoiCells[cellId][i]->getId()];
        if (dist > highestDist) {
            highestDist = dist;
            highestDistListID = i;
        }
    }

    return highestDistListID;
}

#include "blossomgraphmatching.h"

#include <QList>

#include <QDebug>

QVector<int> BlossomGraphMatching::doBlossomGraphMatching(int numberOfNodes, QVector<QPair<int, int> > edges)
{
    const int numberOfEdges = edges.size();

    QVector<QVector<int> > neighbors(numberOfNodes);
    QVector<QVector<int> > incidentEdges(numberOfNodes);

    for (int i = 0; i < numberOfEdges; i++) {
        QPair<int, int> edge = edges[i];
        incidentEdges[edge.first] << i;
        neighbors[edge.first] << edge.second;

        incidentEdges[edge.second] << i;
        neighbors[edge.second] << edge.first;
    }


    //Construct initial matching
    QVector<bool> inMatching(numberOfEdges, false);
    inMatching[0] = true;

    for (int i = 0; i < numberOfEdges; i++) {
        bool neighborEdgeIsInMatching = false;
        foreach (int incidentEdge, incidentEdges[edges[i].first])
            if (inMatching[incidentEdge])
                neighborEdgeIsInMatching = true;

        foreach (int incidentEdge, incidentEdges[edges[i].second])
            if (inMatching[incidentEdge])
                neighborEdgeIsInMatching = true;

        inMatching[i] = !neighborEdgeIsInMatching;
    }

    QVector<int> matching;
    for (int i = 0; i < numberOfEdges; i++) {
        if (inMatching[i])
            matching << i;
    }

    Graph g;
    g.edges = edges;
    for (int i = 0; i < numberOfNodes; i++)
        g.nodes << i;
    g.nodeSubset = QVector<bool>(numberOfNodes, true);

    while (matching.size() < numberOfNodes/2) {
        QVector<int> pathVertices = BlossomGraphMatching::findAugmentingPath(g, inMatching);

        if (pathVertices.size() > 1) {
            //change matching along path
            for (int i = 0; i < pathVertices.size() - 1; i++) {
                int u = pathVertices.at(i);
                int v = pathVertices.at(i+1);

                int edgeId = -1;
                foreach(int incidentEdgeId, incidentEdges[u])
                    if (edges[incidentEdgeId].first == v || edges[incidentEdgeId].second == v)
                        edgeId = incidentEdgeId;

                inMatching[edgeId] = !inMatching[edgeId];
            }

            matching.clear();
            for (int i = 0; i < numberOfEdges; i++) {
                if (inMatching[i])
                    matching << i;
            }

        } else {
            //no augmented path found -> matching is maximal
            qDebug() << "ERROR: maximum matching is not complete!";
            return QVector<int>();
        }

    }

    return matching;
}

QVector<int> BlossomGraphMatching::findAugmentingPath(Graph &g, QVector<bool> &matched)
{
    const int numberOfTotalNodes = g.nodeSubset.size();
    const int numberOfEdges = g.edges.size();
    QVector<Tree> forest;
    QVector<bool> nodeMarked(numberOfTotalNodes, false);
    QVector<bool> edgeMarked(matched);

    QVector<QVector<int> > neighbors(numberOfTotalNodes);
    QVector<QVector<int> > incidentEdges(numberOfTotalNodes);

    for (int i = 0; i < numberOfEdges; i++) {
        QPair<int, int> edge = g.edges[i];
        incidentEdges[edge.first] << i;
        neighbors[edge.first] << edge.second;

        incidentEdges[edge.second] << i;
        neighbors[edge.second] << edge.first;
    }

    for (int i = 0; i < numberOfTotalNodes; i++) {
        if (g.nodeSubset[i]) {
            bool isExposed = true;
            for (int j = 0; j < incidentEdges[i].size(); j++) {
                if (matched[incidentEdges[i][j]]) {
                    isExposed = false;
                }
            }
            if (isExposed) {
                Tree t;
                t.nodeSubset = QVector<bool>(numberOfTotalNodes, false);
                t.nodeSubset[i] = true;
                t.distances = QVector<int>(numberOfTotalNodes, -1);
                t.distances[i] = 0;
                t.predecessors = QVector<int>(numberOfTotalNodes, -1);
                t.predecessors[i] = i;
                t.nodes << i;
                t.root = i;
                forest << t;
            }
        }
    }

    //Find first v
    int v = -1;
    int vTree = -1;
    if (forest.size() > 0) {
        v = forest[0].root;
        vTree = 0;
    } else {
        qDebug() << "ERROR: No unmarked vertex after intialization";
    }

    while (v != -1) {
        //find first unmarked edge containing v
        int edgeId = -1;
        for (int i = 0; i < incidentEdges[v].size() && edgeId == -1; i++)
            if (!edgeMarked[incidentEdges[v][i]])
                edgeId = incidentEdges[v][i];

        while (edgeId != -1) {
            int w = g.edges[edgeId].first;
            if (w == v)
                w = g.edges[edgeId].second;

            //find tree containing w
            int wTree = -1;
            for (int i = 0; i < forest.size() && wTree == -1; i++)
                if (forest[i].nodeSubset[w])
                    wTree = i;

            if (wTree == -1) {
                //w is matched (i.e. not contained in a tree of the forest)
                //find matched edge and matched node
                int matchedEdge = -1;
                for (int i = 0; i < numberOfEdges; i++)
                    if (matched[i] && (g.edges[i].first == w || g.edges[i].second == w))
                        matchedEdge = i;

                int x = g.edges[matchedEdge].first;
                if (x == w)
                    x = g.edges[matchedEdge].second;

                //add w and matched vertex x to tree of v
                forest[vTree].nodes << w;
                forest[vTree].nodeSubset[w] = true;
                forest[vTree].predecessors[w] = v;
                forest[vTree].distances[w] = forest[vTree].distances[v] + 1;
                forest[vTree].edgeIds << edgeId;

                forest[vTree].nodes << x;
                forest[vTree].nodeSubset[x] = true;
                forest[vTree].predecessors[x] = w;
                forest[vTree].distances[x] = forest[vTree].distances[w] + 1;
                forest[vTree].edgeIds << matchedEdge;
            } else {
                //if distance(w, root(w)) is odd, do nothing
                if (forest[wTree].distances[w] % 2 == 0) {
                    if (vTree != wTree) {
                        //found an augmenting patch
                        QList<int> path;
                        int node = v;
                        path << v;
                        while (node != forest[vTree].root) {
                            node = forest[vTree].predecessors[node];
                            path.prepend(node);
                        }
                        node = w;
                        path << node;
                        while (node != forest[wTree].root) {
                            node = forest[wTree].predecessors[node];
                            path << node;
                        }

                        QVector<int> pathVect = path.toVector();
                        return pathVect;
                    } else {
                        //found a blossom (cycle with 2k+1 edges)

//                        //construct path from v to root(v)
//                        QVector<int> predecessorsV;
//                        QVector<int> blossom;
//                        int node = v;
//                        predecessorsV << node;
//                        while (node != forest[vTree].root) {
//                            node = forest[vTree].predecessors[node];
//                            predecessorsV << node;
//                        }

//                        //now, start constructing a path from w to root(w) and abort, once the path intersects with the first path
//                        node = w;
//                        blossom << node;
//                        int predIndex = predecessorsV.indexOf(node);
//                        while (node != forest[vTree].root && predIndex == -1) {
//                            node = forest[vTree].predecessors[node];
//                            blossom << node;
//                            predIndex = predecessorsV.indexOf(node);
//                        }

//                        for (int i = predIndex-1; i > 0; i--) {
//                            blossom << predecessorsV[i];
//                        }

//                        //Construct contracted graph and matching
//                        Graph gContracted;
//                        gContracted.nodeSubset = g.nodeSubset;
//                        foreach (int vertex, blossom) {
//                            gContracted.nodeSubset[vertex] = false;
//                        }
//                        gContracted.nodeSubset << true;
//                        for (int i = 0; i < gContracted.nodeSubset; i++)
//                            if (gContracted.nodeSubset[i])
//                                gContracted.nodes << i;

//                        int newNumberOfTotalNodes = numberOfTotalNodes + 1;
//                        QVector<bool> matchedContracted;
//                        for (int i = 0; i < numberOfEdges; i++) {
//                            QPair<int, int> edge = g.edges[i];
//                            bool firstNodeExists = gContracted.nodeSubset[edge.first];
//                            bool secondNodeExists = gContracted.nodeSubset[edge.second];
//                            if (firstNodeExists && secondNodeExists) {
//                                gContracted.edges << g.edges[i];
//                                matchedContracted << matched[i];
//                            } else if (firstNodeExists) {
//                                gContracted.edges << QPair<int, int>(edge.first, newNumberOfTotalNodes - 1);
//                                matchedContracted << matched[i];
//                            } else if (secondNodeExists) {
//                                gContracted.edges << QPair<int, int>(edge.second, newNumberOfTotalNodes - 1);
//                                matchedContracted << matched[i];
//                            }
//                        }

//                        //find augmenting path in contracted graph
//                        QVector<int> contractedPath = findAugmentingPath(gContracted, matchedContracted);

//                        //TODO lift path to old graph
//                        //TODO return path
                    }
                }
            }

            //mark current edge
            edgeMarked[edgeId] = true;

            //Find new unmarked edge containing v
            edgeId = -1;
            for (int i = 0; i < incidentEdges[v].size() && edgeId == -1; i++)
                if (!edgeMarked[incidentEdges[v][i]])
                    edgeId = incidentEdges[v][i];
        }

        //mark v
        nodeMarked[v] = true;

        //Find new v
        v = -1;
        for (int i = 0; i < forest.size() && v == -1; i++)
            foreach(int potentialV, forest[i].nodes)
                if (forest[i].distances[potentialV] % 2 == 0 && !nodeMarked[potentialV]) {
                    v = potentialV;
                    vTree = i;
                    break;
                }
    }

    return QVector<int>();
}

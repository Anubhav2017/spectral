#ifndef BLOSSOMGRAPHMATCHING_H
#define BLOSSOMGRAPHMATCHING_H

#include <QVector>
#include <QPair>

struct Graph {
    QVector<bool> nodeSubset;
    QVector<int> nodes;
    QVector<QPair<int, int> > edges;
};

struct Tree {
    QVector<bool> nodeSubset;
    QVector<int> distances;
    QVector<int> predecessors;
    QVector<int> nodes;
    QVector<int> edgeIds;
    int root;
};

class BlossomGraphMatching
{
public:
    static QVector<int> doBlossomGraphMatching(int numberOfNodes, QVector<QPair<int, int> > edges);
    static QVector<int> findAugmentingPath(Graph &g, QVector<bool> &matched);
};

#endif // BLOSSOMGRAPHMATCHING_H

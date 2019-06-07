
///
/// @file DataIO.h
/// @author Petra Gospodnetic
///
/// Functions used for position list handling
///

// TODO: Rename DataIO.h to something more appropriate insinuating the use for position lists.

#ifndef DATAIO_H
#define DATAIO_H

#include <vector>
#include <string>

bool loadPositions( std::vector< std::vector<double> >& positions, std::string filename);

bool sortPositions(std::vector< std::vector<double> >& readPosVector, std::vector< std::vector<double> >& topPositions, std::vector< std::vector<double> >& sidePositions);

bool writePositions(std::vector< std::vector<double> >& positionList, std::string filename);

#endif

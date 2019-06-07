
///
/// @file DataIO.cpp
/// @author Petra Gospodnetic
///

#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DataIO.h"

#define BUFF_SIZE 512

// TODO: add the append argument. If the function is called twice on the same vector
//          should the positions be appended to the current positions or the vector
//          should be cleared first?
bool loadPositions( std::vector< std::vector<double> >& positions, std::string filename)
{
    using namespace std;

    std::cout << "Loading values from file" << endl;

    ifstream file(filename.c_str());

    if (!file.is_open())
    {
        cerr << "Error! Could not access the file" << endl;
        return false;
    }

    // parse .txt file
    // skip blank rows - TODO
    // compare first character in row with #
    // # represents row which defines Position number
    // check if the first character in the following row is a number - TODO
    // mean value vector follows in the next row

    {
        char check_line = 'X';

        char line[BUFF_SIZE];
        char meanLine[BUFF_SIZE];

        while (file.get(check_line))
        {

            file.getline(line, BUFF_SIZE);
            file.getline(meanLine, BUFF_SIZE);

            // split the line into separate numbers
            char* tokens;

            tokens = strtok(meanLine, " ");
            vector <double> tmpPositions;
            while (tokens != NULL)
            {
                tmpPositions.push_back(atof(tokens));
                tokens = strtok(NULL, " ");
            }
            positions.push_back(tmpPositions);
        }
    }
    cout << "Number of positions: " << positions.size() << endl;

    return true;
}

bool writePositions(std::vector< std::vector<double> >& positionList, std::string filename)
{
    using namespace std;
    ofstream file;
    file.open(filename.c_str());

    if (!file.is_open())
    {
        cerr << "Error! File could not be opened!" << endl;
        return false;
    }

    for (unsigned int i = 0; i < positionList.size(); i++)
    {
        file << "#Position " << i << ": " << endl;
        for(unsigned int j = 0; j < positionList[0].size(); j++ )
        {
            file << positionList[i][j] << " ";
        }

        file << endl;
    }

    file.close();
    return true;
}

bool sortPositions(std::vector< std::vector<double> >& readPosVector, std::vector< std::vector<double> >& topPositions, std::vector< std::vector<double> >& sidePositions)
{
    for (unsigned int i = 0; i < readPosVector.size(); i++)
    {
        double dir_x = fabs(round(readPosVector[i][3]));
        double dir_y = fabs(round(readPosVector[i][4]));
        double dir_z = fabs(round(readPosVector[i][5]));


        if(dir_x == 0.0 && dir_y == 0.0 && dir_z == 1.0)
        {
            // Camera looks downwards
            topPositions.push_back(readPosVector[i]);
        }
        else
        {
            sidePositions.push_back(readPosVector[i]);
        }
    }
    return true;
}

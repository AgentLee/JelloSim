#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

#include "particles.h"

class Tetrahedron
{
public:
	Tetrahedron();

    int numTetra;
    std::vector<Eigen::Matrix<uint, 4, 1>> particleIndices;

    // Reads Tetgen file and stores data
    void tetgen_readLine(std::ifstream &fin, int nodesPerTet);
	void tetgen_readEleFile(const std::string &inputFileName);
};

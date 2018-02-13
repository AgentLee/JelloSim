#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

class Tetrahedron
{
public:
	Tetrahedron();

    int numTetra;
    std::vector<Eigen::Matrix<uint, 4, 1>> particleIndices;

	void readEle(const std::string &inputFileName);
};

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include <stdlib.h>
#include <time.h>

#include "particles.h"

using T = double;

class Tetrahedrons
{
public:
	Tetrahedrons();

    int numTetra;
    std::vector<Eigen::Matrix<uint, 4, 1>> particleIndices;

    std::vector<T> undeformedVolume; //W
    std::vector<Eigen::Matrix<T, 3, 3>> restDeformation; //Dm
    std::vector<Eigen::Matrix<T, 3, 3>> restInverseDeformation; //Dm inverse

    void computeRestDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices ); // Calculate rest Deformation
    void computeInvRestDeformation( uint tetraIndex ); // Calculate inverse rest Deformation
    void computeUndeformedVolume( uint tetraIndex ); // Calculate undeformed Volume

    void computeNewDeformation( Eigen::Matrix<T, 3, 3> newDef, uint tetraIndex ); // Calculate new deformation

    // Reads Tetgen file and stores data
    void tetgen_readLine(std::ifstream &fin, int nodesPerTet);
	void tetgen_readEleFile(const std::string &inputFileName);
};

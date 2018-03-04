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
	Tetrahedrons( float density, float _poissonsRatio, float Y );

    int numTetra;
    float density;
    float poissonsRatio;
    float youngsModulus;
    std::vector<Eigen::Matrix<uint, 4, 1>> particleIndices;

    //Needed for FEM Calculations
    std::vector<T> undeformedVolume; //W
    std::vector<Eigen::Matrix<T, 3, 3>> restDeformation; //Dm
    std::vector<Eigen::Matrix<T, 3, 3>> restInverseDeformation; //Dm inverse
    std::vector<Eigen::Matrix<T, 3, 3>> undefVol_into_restInvDefTranspose; //undeformedVolume[tetraIndex] * restInverseDeformation[tetraIndex].transpose()

    // Precompute
    void computeRestDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices ); // Calculate rest Deformation
    void computeInvRestDeformation( uint tetraIndex ); // Calculate inverse rest Deformation
    void computeUndeformedVolume( uint tetraIndex ); // Calculate undeformed Volume
    void computeUndefVol_into_restInvDefTranspose( uint tetraIndex ); // Calculate undeformed Volume( uint tetraIndex ); // Calculate undeformed Volume

    void addMass( uint tetraIndex, std::shared_ptr<Particles> vertices );

    // Main Update Loop
    Eigen::Matrix<T, 3, 3> computeNewDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices ); // Calculate new deformation
    Eigen::Matrix<T, 3, 3> computeF( uint tetraIndex, Eigen::Matrix<T, 3, 3>& Ds ); // Calculate F
    Eigen::Matrix<T, 3, 3> computeP( uint tetraIndex, const Eigen::Matrix<T, 3, 3>& F, int frame, bool &collided ); // Calculate P
    Eigen::Matrix<T, 3, 3> computeH( uint tetraIndex, Eigen::Matrix<T, 3, 3>& P ); // Calculate H

    //Add a force to every particle that comprises the tetrahedron
    void addForces( uint tetraIndex, std::shared_ptr<Particles> vertices, Eigen::Matrix<T,3,3>& H );

    // Reads Tetgen file and stores data
    void tetgen_readLine(std::ifstream &fin, int nodesPerTet);
	void tetgen_readEleFile(const std::string &inputFileName);
};

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

#include "bounds.h"

using T = double;
constexpr int dim = 3;
const int num = 50;
const double deltaT = 0.1;

class Particles
{
public:
	int numParticles;
	std::vector<T> mass;
	std::vector<Eigen::Matrix<T, 3, 1>> vel;
	std::vector<Eigen::Matrix<T, 3, 1>> pos;
	std::vector<Eigen::Matrix<T, 3, 1>> force;

	Particles(int n, T initialMass);
	void updateAllParticlePositions(T dt);
	void updateAllParticleVelocities(T dt, Bounds& FixedRegion);
	void updateParticlePosition(T dt, uint index);
	void updateParticleVelocity(T dt, uint index);

	// Reads Tetgen file and stores data
	void tetgen_readLine(std::ifstream &fin, int numDims, int numAtt);
	void tetgen_readNode(const std::string &inputFileName);

	void writePartio(std::string particleFile);
};

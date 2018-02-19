#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

using T = double;

class Particles
{
public:
	int numParticles;
	std::vector<T> mass;
	std::vector<Eigen::Matrix<T,3,1>> vel;
	std::vector<Eigen::Matrix<T,3,1>> pos;
	std::vector<Eigen::Matrix<T,3,1>> force;

	Particles(int n, T initialMass);
	void initializeParticles_random();
	void updateParticlePositions(T dt);
	void updateParticleVelocity(T dt);

	// Reads Tetgen file and stores data
	void tetgen_readLine(std::ifstream &fin, int numDims, int numAtt);
	void tetgen_readNode(const std::string &inputFileName);
};

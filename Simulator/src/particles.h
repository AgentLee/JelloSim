#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>

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

public:
	Particles(int n, T initialMass);
	void initializeParticles_random();
	void updateParticlePositions(T dt);
};
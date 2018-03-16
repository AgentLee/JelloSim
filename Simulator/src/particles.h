#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

#include <GL/glew.h>
#include <GLM/glm.hpp>
#include <GLM/gtc/matrix_transform.hpp>
#include <GLM/gtc/type_ptr.hpp>

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
	void updateAllParticleVelocities(T dt);
	void updateParticlePosition(T dt, uint index);
	void updateParticleVelocity(T dt, uint index);

	void updateAllParticlePositionsBackwardEuler(T dt);
	void updateAllParticleVelocitiesBackwardEuler(T dt);
	void updateParticlePositionBackwardEuler(T dt, uint index);
	void updateParticleVelocityBackwardEuler(T dt, uint index);

	// Reads Tetgen file and stores data
	void tetgen_readLine(std::ifstream &fin, int numDims, int numAtt);
	void tetgen_readNode(const std::string &inputFileName);

	void writePartio(std::string particleFile);

	void computeForce(int index);

	const T Ks = 1.0;
	const T Kd = 1.0;

	T mass_damping = 2.0;
	glm::vec3 X, V, F, Vnew;
	glm::mat3 kMat, mMat; 
	std::vector<glm::mat3> dForceDX;      // Implicit Modified Backward Euler, ∂F/∂X, Partial Differential (ΔForces / ΔPositions)
	std::vector<glm::mat3> dForceDV;      // Implicit Modified Backward Euler, ∂F/∂V, Partial Differential (ΔForces / ΔVelocitys)
};

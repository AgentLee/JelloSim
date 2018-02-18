#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <time.h>

#include <memory>
#include <unordered_map>

#include "tetrahedrons.h"
#include "particles.h"
#include "sdf.h"

using T = double;

class Sim
{
public:
	std::shared_ptr<Tetrahedrons> tetras;
	std::shared_ptr<Particles> vertices;
	
	//hash tables for collision forces and for all other forces --> 2 tables so as to split it up 
	//into 'forces on some of the particles' versus 'forces on every particle'
	std::unordered_map<long, Eigen::Matrix<T,3,1> > collisionForceMap;
	std::unordered_map<long, Eigen::Matrix<T,3,1> > elasticForceMap;

	Sim(std::shared_ptr<Tetrahedrons> tetrahedronList, std::shared_ptr<Particles> particleList);

	void init();
	void clean();
	bool checkHashElasticForces( int i, int j, Eigen::Matrix<T,3,1>& force );
	void computeElasticForce( int tetraIndex );
	void eulerIntegration();
	void update();
	void checkCollisions();
};
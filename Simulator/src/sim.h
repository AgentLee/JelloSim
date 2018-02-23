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

	Sim(std::shared_ptr<Tetrahedrons> tetrahedronList, std::shared_ptr<Particles> particleList);

	void init();
	void clean();
	void computeElasticForces( int tetraIndex, int frame );
	void eulerIntegration(float dt);
	void update(float dt, int i);
	void checkCollisions(float dt);
};
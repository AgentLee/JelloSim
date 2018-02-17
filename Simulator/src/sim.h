#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <time.h>

#include <memory>

#include "tetrahedron.h"
#include "particles.h"

class Sim
{
public:
	std::shared_ptr<Tetrahedron> tetras;
	std::shared_ptr<Particles> vertices;

	Sim(std::shared_ptr<Tetrahedron> tetrahedronList, std::shared_ptr<Particles> particleList);

	void init();
	void eulerIntegration();
	void update();
	void checkCollisions();
};
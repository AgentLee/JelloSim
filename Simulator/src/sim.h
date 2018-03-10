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
#include "mesh.h"
#include "sdf.h"
#include "triangles.h"

using T = double;

class Sim
{
public:
	std::vector<std::shared_ptr<Mesh>> MeshList;

	Sim( std::vector<std::shared_ptr<Mesh>>& MeshList );

	void init();
	void clean();
	
	void eulerIntegration(float dt);
	void checkCollisions(float dt, bool &collided);

	void addExternalForces();
	void computeElasticForces( int frame, bool &collided );

	void fixParticlePosition( Eigen::Matrix<T, 3, 1>& particleVel, Eigen::Matrix<T, 3, 1>& particlePos );

	void update(float dt, int frame, bool &collided);

	int colliding(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, float *t);
};
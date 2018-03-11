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
#include "bounds.h"
#include "sdf.h"

using T = double;

class Sim
{
public:
	std::vector<std::shared_ptr<Mesh>> MeshList;

	Sim( std::vector<std::shared_ptr<Mesh>>& MeshList );

	void init();
	void clean();
	
	void eulerIntegration(float dt);
	void eulerIntegrationWithCollisionTesting(float dt, bool& collided);
	void checkCollisions(float dt, uint j, bool& collided);

	void addExternalForces();
	void computeElasticForces( int frame, bool &collided );

	void fixParticlePosition( Eigen::Matrix<T, 3, 1>& particleVel, Eigen::Matrix<T, 3, 1>& particlePos );
	void resolveCollisions( std::shared_ptr<Triangles>& triangles, std::shared_ptr<Particles>& vertices, 
							Intersection& isect, Eigen::Matrix<T, 3, 1>& displacement,
							Eigen::Matrix<T, 3, 1>& particlePos, Eigen::Matrix<T, 3, 1>& particleVel );

	void update(float dt, int frame, bool &collided);
};
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
#include "triangles.h"

using T = double;

class Sim
{
public:
	std::vector<std::shared_ptr<Mesh>> MeshList;

	Sim( std::vector<std::shared_ptr<Mesh>>& MeshList );

	void init();
	void clean();

	void reComputeMeshAttributes();
	
	void eulerIntegration(float dt);
	void eulerIntegrationWithCollisionTesting(float dt);

	void addExternalForces();
	void computeElasticForces( int frame );

	void SDF_Collisions(float dt, uint j);
	void fixParticlePosition( Eigen::Matrix<T, 3, 1>& particleVel, Eigen::Matrix<T, 3, 1>& particlePos );
	void Mesh_Collisions(float dt, uint i, uint j, std::shared_ptr<Particles>& vertices, int particleIndex, bool& collided);
	void resolveCollisions( std::shared_ptr<Triangles>& triangles, std::shared_ptr<Particles>& vertices,
							Intersection& isect, Eigen::Matrix<T, 3, 1>& displacement,
							Eigen::Matrix<T, 3, 1>& particlePos, Eigen::Matrix<T, 3, 1>& particleVel );

	void update(float dt, int frame);

    bool LineTriangleIntersection(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, Intersection *isect);
};
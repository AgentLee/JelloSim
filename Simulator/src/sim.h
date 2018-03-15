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
    Bounds FixedRegion;

	Sim( std::vector<std::shared_ptr<Mesh>>& MeshList, Bounds& FixedRegion );

	void init();
	void clean();

	void reComputeMeshAttributes();
	void computeAllForces( int frame );

    void eulerIntegration(float dt);
	void eulerIntegrationWithSDF_Collisions(float dt);

	void update(float dt, int frame);

	void SDF_Collisions(float dt, uint j);
    void Collisions(float dt);
	
	bool IntersectionTesting(std::shared_ptr<Mesh> meshA, std::shared_ptr<Mesh> meshB, int vertexID, Intersection *isect);
};
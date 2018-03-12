/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SET_POSITIONS 1
#define SET_VELOCITIES 1
#define PAPER 0
#define INTER_OBJECT_COLLISIONS 1

Sim::Sim( std::vector<std::shared_ptr<Mesh>>& MeshList ) : MeshList(MeshList)
{}

void Sim::init()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->initMeshForSim();
	}
}

void Sim::clean()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->clearForces();
	}
}

void Sim::reComputeMeshAttributes()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->triangles->computeNormals(MeshList[i]->vertices);
		MeshList[i]->calcBounds();
	}
}

void Sim::computeAllForces( int frame )
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		//////////------------------ EXTERNAL FORCES ---------------/////////////
		//Not abstracting to mesh because we could add forces that only operate on 
		//some of the vertices and things like that which is annoying to genarlize for
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		
		for(int j=0; j<vertices->numParticles; j++)
		{
			vertices->force[j](1) -= 9.81 * vertices->mass[j]; // gravity
		}

		//////////------------------ ELASTIC FORCES ---------------/////////////
		MeshList[i]->computeElasticForcesOnMesh(frame);
	}
}

void Sim::update(float dt, int frame)
{
	clean(); //clears forces
    reComputeMeshAttributes(); //compute triangle normals for all meshes
	computeAllForces(frame); //computes and adds elastic forces to each particle

#if INTER_OBJECT_COLLISIONS
    //eulerIntegrationWithCollisionTesting(dt);
    eulerIntegration(dt);
    Collisions(dt);
#else
    eulerIntegrationWithSDF_Collisions(dt);
#endif
}

void Sim::eulerIntegration(float dt)
{
    for(uint i=0; i<MeshList.size(); i++)
    {
        MeshList[i]->setPrevVerts();
        std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
        vertices->updateAllParticleVelocities(dt);
        vertices->updateAllParticlePositions(dt);
    }
}

void Sim::eulerIntegrationWithSDF_Collisions(float dt)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
        SDF_Collisions(dt, i);
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		vertices->updateAllParticleVelocities(dt);
		vertices->updateAllParticlePositions(dt);
	}
}

// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
bool Sim::LineTriangleIntersection(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, Intersection *isect)
{
	// create ray and call triangle intersection for all triangles..
	// store triangle reference
	Ray r;
	r.origin = origPos;
	r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);
    float length = r.direction.norm();
    r.direction.normalize();

	// assuming only 2 meshes for now..
	int tris = MeshList[1]->triangles->triFaceList.size();
	isect->t = std::numeric_limits<T>::infinity();
    isect->hit = false;
	isect->triangleIndex = -1;
	for(int i=0; i < tris; i++)
	{
		float tTemp = isect->t;
        Eigen::Matrix<T, 3, 1> baryCoords;
		bool intersect = MeshList[1]->triangles->intersect(r, i, &tTemp, MeshList[1]->vertices, &baryCoords);
		if(intersect && isect->t > tTemp && tTemp <= length)
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp;
            isect->BarycentricWeights = baryCoords;
            isect->normal = MeshList[1]->triangles->triNormalList[i];
            isect->point = r.origin + tTemp * r.direction;
		}
	}

	return isect->hit;
}

void Sim::SDF_Collisions(float dt, uint j)
{
	// Check if the Mesh hit the ground or any other solid piece of geometry that is 
	// essentially an infinite Mass Rigid Body that doesnt move
	std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
	std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

	for(int i = 0; i< vertices->numParticles; ++i)
	{
		Eigen::Matrix<T, 3, 1> p(vertices->pos[i]);

		// Transform the vertex to the plane's local space
		// Assume the plane is at the origin
		Eigen::Matrix<T, 4, 1> n = Eigen::Matrix<T, 4, 1>(0, 1, 0, 0);
		float sdf = SDF::sdPlane(p, n);

		// Check if particle went through the surface
		if(sdf < 0) 
		{
			if( vertices->vel[i].dot(Eigen::Matrix<T, 3, 1>(0, 1, 0)) < 0 )
			{
				vertices->vel[i] = Eigen::Matrix<T, 3, 1>::Zero();
			}
		}
	}
}

void Sim::Collisions(float dt)
{
    for(uint i=0; i<MeshList.size(); i++)
    {
        SDF_Collisions(dt, i);
    }

    std::vector<std::vector<Eigen::Matrix<T, 3, 1>>> newPositions;
    std::vector<std::vector<Eigen::Matrix<T, 3, 1>>> newVelocities;
    for(uint i = 0; i < MeshList.size(); i++)
    {
        std::vector<Eigen::Matrix<T, 3, 1>> pos;
        std::vector<Eigen::Matrix<T, 3, 1>> vel;
        for(int j = 0; j < MeshList[i]->vertices->numParticles; j++)
        {
            pos.push_back(MeshList[i]->vertices->pos[j]);
            vel.push_back(MeshList[i]->vertices->vel[j]);
        }
        newPositions.push_back(pos);
        newVelocities.push_back(vel);
    }


    for(uint i=0; i<MeshList.size(); i++)
    {
        for (uint j = 0; j < MeshList.size(); j++)
        {
            std::shared_ptr<Mesh> fallingMesh = MeshList[i];
            std::shared_ptr<Mesh> standingMesh = MeshList[j];

            bool intersects = Intersect_AABB_with_AABB(fallingMesh->AABB, standingMesh->AABB);

            if (intersects) {
                for (int k = 0; k < fallingMesh->vertices->numParticles; ++k) {
                    std::shared_ptr<Triangles> triangles = standingMesh->triangles; //i --> current mesh's triangles

                    // for every vertex create a ray from the current position to its projected position in the next frame
                    // See if that ray intersects any triangle that Belongs to the other mesh.
                    //Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[particleIndex] + vertices->vel[particleIndex] * dt;
                    Intersection isect;
                    LineTriangleIntersection(fallingMesh->vertices->pos[k], fallingMesh->prevVerts->pos[k], &isect);

                    if (isect.hit) {
                        //resolve collisions between vertex and a triangle
                        Eigen::Matrix<T, 3, 1> displacement;
                        displacement = fallingMesh->vertices->pos[k] - isect.point;

                        // set position or velocity or both OR use momentum OR use paper's implementation with normal reaction forces and friction
                        Eigen::Matrix<uint, 3, 1> verticesOfTriangle;
                        verticesOfTriangle = triangles->triFaceList[isect.triangleIndex];

                        if (SET_POSITIONS) {
                            //---------------- Setting Positions ---------------------------------
                            //fallingMesh->vertices->pos[k] = fallingMesh->prevVerts->pos[k]; //isect.point.cast<T>() - displacement; // 1/4th to moving vertex
                            newPositions[i][k] = isect.point - 1.3 * displacement;
                            //displacement = 0.75f*displacement;				// 3/4th to moving triangle

                            // move vertices of triangle according to barycentric weights
//                            standingMesh->vertices->pos[verticesOfTriangle[0]] +=
//                                    isect.BarycentricWeights[0] * displacement;
//                            standingMesh->vertices->pos[verticesOfTriangle[1]] +=
//                                    isect.BarycentricWeights[1] * displacement;
//                            standingMesh->vertices->pos[verticesOfTriangle[2]] +=
//                                    isect.BarycentricWeights[2] * displacement;

                            //        std::cout<< "disp: " << displacement << std::endl;
                            //        std::cout<< "0: " << isect.BarycentricWeights[0] << std::endl;
                            //        std::cout<< "1: " << isect.BarycentricWeights[1] << std::endl;
                            //        std::cout<< "2: " << isect.BarycentricWeights[2] << std::endl;

                        }
                        if (SET_VELOCITIES) {
                            //fallingMesh->vertices->vel[k] = Eigen::Matrix<T, 3, 1>::Zero();
                            newVelocities[i][k] = Eigen::Matrix<T, 3, 1>::Zero();
                            //MeshList[0]->vertices->force[i](1) += 9.81 * MeshList[1]->vertices->mass[i]; // gravity
//
//                            MeshList[1]->vertices->vel[verticesOfTriangle[0]] = Eigen::Matrix<T, 3, 1>::Zero();
//                            MeshList[1]->vertices->vel[verticesOfTriangle[1]] = Eigen::Matrix<T, 3, 1>::Zero();
//                            MeshList[1]->vertices->vel[verticesOfTriangle[2]] = Eigen::Matrix<T, 3, 1>::Zero();
                        }
                    }
                }
            }
        }
    }

    for(uint i = 0; i < MeshList.size(); i++)
    {
        for(int j = 0; j < MeshList[i]->vertices->numParticles; j++)
        {
            MeshList[i]->vertices->pos[j] = newPositions[i][j];
            MeshList[i]->vertices->vel[j] = newVelocities[i][j];
        }
    }
}

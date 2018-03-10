/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SPRING_COLLISION 0
#define SET_POSITIONS 1
#define SET_VELOCITIES 0
#define PAPER 0

Sim::Sim( std::vector<std::shared_ptr<Mesh>>& MeshList ) : MeshList(MeshList)
{}

void Sim::init()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

		// Precompute rest deformation (Dm), volume, inverse Dm, and volume*inverseDmTranspose for each tetrahedron
		tetras->restDeformation.resize(tetras->numTetra);
		tetras->restInverseDeformation.resize(tetras->numTetra);
		tetras->undeformedVolume.resize(tetras->numTetra);
		tetras->undefVol_into_restInvDefTranspose.resize(tetras->numTetra);

		for(int i=0; i<tetras->numTetra; i++)
		{
			tetras->computeRestDeformation( i, vertices );
			tetras->computeInvRestDeformation( i );
			tetras->computeUndeformedVolume( i );
			tetras->computeUndefVol_into_restInvDefTranspose( i );

			tetras->addMass( i, vertices );
		}
	}
}

void Sim::clean()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
		//Set forces for all vertices/particles to zero
		for( int i=0; i < vertices->numParticles; i++ )
		{
			vertices->force[i] = Eigen::Matrix<T, 3, 1>::Zero();
		}
	}
}

void Sim::eulerIntegration(float dt)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;

		vertices->updateParticleVelocity(dt);
		vertices->updateParticlePositions(dt);
	}
}

void Sim::eulerIntegrationWithCollisionTesting(float dt, bool& collided)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		checkCollisions(dt, i, collided); //Apply Forces to particles that occur through collision
		
		if(collided)
		{
			vertices->updateParticleVelocity(dt);
			vertices->updateParticlePositions(dt);
		}
	}
}

void Sim::addExternalForces()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
		
		for(int i=0; i<vertices->numParticles; i++)
		{
			vertices->force[i](1) -= 9.81 * vertices->mass[i]; // gravity
		}
	}
}

void Sim::computeElasticForces( int frame, bool& collided )
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

		// Loop through tetras
		for(int tetraIndex=0; tetraIndex < tetras->numTetra; tetraIndex++)
		{
			Eigen::Matrix<T,3,3> newDeformation = tetras->computeNewDeformation( tetraIndex, vertices ); // Compute Ds, the new deformation
			Eigen::Matrix<T,3,3> F 				= tetras->computeF( tetraIndex, newDeformation ); // Compute F = Ds(Dm_inv)
			Eigen::Matrix<T,3,3> P 				= tetras->computeP( tetraIndex, F, frame, collided ); // Compute Piola (P)
			Eigen::Matrix<T,3,3> H 				= tetras->computeH( tetraIndex, P ); // Compute Energy (H)

			tetras->addForces( tetraIndex, vertices, H );// Add energy to forces (f += h)
		}
	}
}

void Sim::update(float dt, int frame, bool& collided)
{
	clean(); //clears forces

    addExternalForces();
	computeElasticForces(frame, collided); //computes and adds elastic forces to each particle

    eulerIntegrationWithCollisionTesting(dt, collided);
}


// TODO: Re-compute normals of faces

// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
bool Sim::LineTriangleIntersection(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, Intersection *isect)
{
	// create ray and call triangle intersection for all triangles..
	// store triangle reference

	Ray r;
	r.origin = origPos;
	r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);

	// assuming only 2 meshes for now..
	int tris = MeshList[1]->triangles->triFaceList.size();
	isect->t = -1;
	isect->triangleIndex = -1;
	for(int i=0; i < tris; i++)
	{
		float tTemp = isect->t;
		bool intersect = MeshList[1]->triangles->intersect(r, i, &tTemp, MeshList[1]->vertices);
		if(intersect && isect->t > tTemp)
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp;
		}
	}

	return isect->hit;
}



void Sim::fixParticlePosition(Eigen::Matrix<T, 3, 1>& particleVel, Eigen::Matrix<T, 3, 1>& particlePos)
{
	if(SPRING_COLLISION)
	{
		// // Apply zero length spring
		// Eigen::Matrix<T, 3, 1> vSurf = p;
		// p[1] -= sdf;

		// // WHERE IS k DEFINED
		// T k = (T)0;
		// vertices->force.at(i) += (-k * (p - vSurf));
	}
	else
	{
		// Move the particle up to the surface
		// Subtract the y component by the SDF	

		// vertices->pos[i](0, 1) = 0.00001;

		if(particleVel.dot(Eigen::Matrix<T, 3, 1>(0, 1, 0)) < 0)
		{
			particleVel = Eigen::Matrix<T, 3, 1>::Zero();
		}
	}
}

void Sim::resolveCollisions( std::shared_ptr<Triangles>& triangles, std::shared_ptr<Particles>& vertices, 
	Intersection& isect, Eigen::Matrix<T, 3, 1>& displacement,
	Eigen::Matrix<T, 3, 1>& particlePos, Eigen::Matrix<T, 3, 1>& particleVel )
{
	// NOTE: triangles and vertices corresspond to the Triangles and Vertices of the other Mesh

	// set position or velocity or both OR use momentum OR use paper's implementation with normal reaction forces and friction
	Eigen::Matrix<uint, 3, 1> verticesOfTriangle = triangles->triFaceList[isect.triangleIndex];

	if(SET_POSITIONS)
	{
		//---------------- Setting Positions ---------------------------------
		particlePos = isect.point.cast<T>() - 0.25f*displacement; // 1/4th to moving vertex
		displacement = 0.75f*displacement;				// 3/4th to moving triangle

		// move vertices of triangle according to barycentric weights
		vertices->pos[verticesOfTriangle[0]] += isect.BarycentricWeights[0] * displacement;
		vertices->pos[verticesOfTriangle[1]] += isect.BarycentricWeights[1] * displacement;
		vertices->pos[verticesOfTriangle[2]] += isect.BarycentricWeights[2] * displacement;
	}
	if(SET_VELOCITIES)
	{
		particleVel = Eigen::Matrix<T, 3, 1>::Zero();

		vertices->vel[verticesOfTriangle[0]] = Eigen::Matrix<T, 3, 1>::Zero();
		vertices->vel[verticesOfTriangle[1]] = Eigen::Matrix<T, 3, 1>::Zero();
		vertices->vel[verticesOfTriangle[2]] = Eigen::Matrix<T, 3, 1>::Zero();
	}
	if(PAPER)
	{
	}
}

void Sim::checkCollisions(float dt, uint j, bool& collided)
{
	//-------- First ---------------
	// Check if any of the Meshes hit the ground or any other solid piece of geometry that is 
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
			collided = true;
			fixParticlePosition( vertices->vel[i], vertices->pos[i] );
		}
	}

	//-------- Second ---------------
	//Check if this mesh's AABB is intersecting with any other mesh's AABB
	for(uint i=0; i<MeshList.size() && i!=j; i++)
	{
		bool intersects = Intersect_AABB_with_AABB( MeshList[j]->AABB, MeshList[i]->AABB );
		
		if(intersects)
		{
			std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
			// Then do a brute force check, i.e loop through all vertice of one mesh 
			// OR use a grid structure or cull the triangles somehow
			
			for(int k = 0; k< vertices->numParticles; ++k)
			{
				std::shared_ptr<Triangles> triangles = MeshList[j]->triangles;

				// for every vertex create a ray from the current position to its projected position in the next frame
				// See if that ray intersects any triangle that Belongs to the other mesh.

				Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[k] + vertices->vel[k] * dt;
				Intersection isect;// = LineTriangleIntersection(triangles, vertices->pos[k], projectedPos);

				if(isect.hit)
				{
					//resolve collisions between vertex and a triangle
					Eigen::Matrix<T, 3, 1> displscement = projectedPos - vertices->pos[i];

					resolveCollisions( MeshList[i]->triangles, MeshList[i]->vertices, isect, displscement, 
									vertices->pos[k], vertices->vel[k] );
					collided = true;
				}
			}
		}
	}
}
/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SPRING_COLLISION 0

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

void Sim::computeElasticForces(int frame, bool &collided)
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

void Sim::update(float dt, int frame, bool &collided)
{
	clean(); //clears forces
    checkCollisions(dt, collided); //Apply Forces to particles that occur through collision

    addExternalForces();
	computeElasticForces(frame, collided); //computes and adds elastic forces to each particle

    eulerIntegration(dt);
}

/*
 * loop over all the triangles
 * 	if normal is facing away,
 * 		skip
 * 	else
 * 		find intersection
 *
 *	sort and get the closest intersection
 */

// TODO: Re-compute normals of faces
// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
int Sim::colliding(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, float *t)
{
	// create ray and call triangle intersection for all triangles..
	// store triangle reference

	Ray r;
	r.origin = origPos;
	r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);

	// assuming only 2 meshes for now..
	int tris = MeshList[1]->triangles->triFaceList.size();
	*t = -1;
	int closestTri = -1;
	for(int i=0; i < tris; i++)
	{
		float tTemp = *t;
		bool intersect = MeshList[1]->triangles->intersect(r, i, &tTemp, MeshList[1]->vertices);
		if(intersect && *t > tTemp)
		{
			closestTri = i;
			*t = tTemp;
		}
	}

	return closestTri;
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

void Sim::resolveCollisions( std::shared_ptr<Triangles>& triangles, Intersection& isect, Eigen::Matrix<T, 3, 1>& particlePos )
{
	// move particle outside triangle
	// or set its velocity in accordance to momentum equations
	// or both of the above


}

void Sim::checkCollisions(float dt, bool &collided)
{
	//-------- First ---------------
	// Check if any of the Meshes hit the ground or any other solid piece of geometry that is 
	// essentially an infinite Mass Rigid Body that doesnt move
	for(uint j=0; j<MeshList.size(); j++)
	{
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
	}

	//-------- Second ---------------
	for(uint j=0; j<MeshList.size(); j++)
	{
		//Check if this mesh's AABB is intersecting with any other mesh's AABB
		for(uint i=0; i<MeshList.size() && i!=j; i++)
		{
			bool intersects = Intersect_AABB_with_AABB( MeshList[j]->AABB, MeshList[i]->AABB );
			
			if(intersects)
			{
				std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
				// Then do a brute force check, i.e loop through all vertice of one mesh 
				// OR use a grid structure or cull the triangles somehow
				
				for(int i = 0; i< vertices->numParticles; ++i)
				{
					std::shared_ptr<Triangles> triangles = MeshList[j]->triangles;

					// for every vertex create a ray from the current position to its projected position in the next frame
					// See if that ray intersects any triangle that Belongs to the other mesh.

					//Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[i] + vertices->vel[i] * dt;
					Intersection isect;// = LineTriangleIntersection(triangles, vertices->pos[i], projectedPos);

					if(isect.hit)
					{
						//resolve collisions between vertex and a triangle
						resolveCollisions( triangles, isect, vertices->pos[i] );
					}
				}
			}
		}
	}
}
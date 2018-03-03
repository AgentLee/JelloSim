/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SPRING_COLLISION false

Sim::Sim(std::shared_ptr<Tetrahedrons> tetrahedronList, std::shared_ptr<Particles> particleList) : tetras(tetrahedronList), vertices(particleList)
{}

void Sim::init()
{
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

		tetras->addMass( i, 1000.0f, vertices );
	}
}

void Sim::clean()
{
	//Set forces for all vertices/particles to zero
	for(int i=0; i<vertices->numParticles; i++)
	{
		vertices->force[i] = Eigen::Matrix<T, 1, 3>::Zero();
	}
}

void Sim::eulerIntegration(float dt)
{        
	vertices->updateParticleVelocity(dt);
    vertices->updateParticlePositions(dt);
}

void Sim::computeElasticForces( int tetraIndex, int frame, bool &collided )
{
	Eigen::Matrix<T,3,3> newDeformation = tetras->computeNewDeformation( tetraIndex, vertices ); // Compute Ds, the new deformation
	Eigen::Matrix<T,3,3> F 				= tetras->computeF( tetraIndex, newDeformation ); // Compute F = Ds(Dm_inv)
	Eigen::Matrix<T,3,3> P 				= tetras->computeP( tetraIndex, F, frame, collided ); // Compute Piola (P)
	Eigen::Matrix<T,3,3> H 				= tetras->computeH( tetraIndex, P ); // Compute Energy (H)

	tetras->addForces( tetraIndex, vertices, H );// Add energy to forces (f += h)
}

void Sim::update(float dt, int frame, bool &collided)
{
	clean(); //clears forces
    checkCollisions(dt, collided); //Apply Forces to particles that occur through collision

    // External forces like gravity
    for(int i=0; i<vertices->numParticles; i++)
    {
        vertices->force[i](0,1) -= 9.81 * vertices->mass[i]; //computes and adds elastic forces to each particle
    }

	// Loop through tetras
	for(int i=0; i<tetras->numTetra; i++)
	{
		computeElasticForces(i, frame, collided); //computes and adds elastic forces to each particle
	}

    eulerIntegration(dt);
}

void Sim::checkCollisions(float dt, bool &collided)
{
	// First do brute force SDF for primitives
	// https://gamedev.stackexchange.com/questions/66636/what-are-distance-fields-and-how-are-they-applicable-to-collision-detection
	// Transform the particles to the static mesh's local space
	// If negative, particle went through the surface.
	// If positive, particle didn't hit.
	// If zero, particle on the surface.

	// Later do interobject collisions with a grid+bounding box or BVH
	for(int i = 0; i< vertices->numParticles; ++i)
	{
		Eigen::Matrix<T, 1, 3> p(vertices->pos[i]);

		// Transform the vertex to the plane's local space
		// Assume the plane is at the origin
		Eigen::Matrix<T, 1, 4> n = Eigen::Matrix<T, 1, 4>(0,1,0,0);
		float sdf = SDF::sdPlane(p, n);

		// Particle went through the surface
		if(sdf < 0) {
			collided = true;

			if(SPRING_COLLISION) {
				// // Apply zero length spring
				// Eigen::Matrix<T, 3, 1> vSurf = p;
				// p[1] -= sdf;

				// // WHERE IS k DEFINED
				// T k = (T)0;
				// vertices->force.at(i) += (-k * (p - vSurf));
			}
			else {
				// Move the particle up to the surface
				// Subtract the y component by the SDF	
				
				// vertices->pos[i](0, 1) = 0.00001;

				if(vertices->vel[i].dot(Eigen::Matrix<T, 1, 3>(0, 1, 0)) < 0) {
					vertices->vel[i] = Eigen::Matrix<T, 1, 3>::Zero();
				}
			}
		}
	}
}
/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

Sim::Sim(std::shared_ptr<Tetrahedrons> tetrahedronList, std::shared_ptr<Particles> particleList) : tetras(tetrahedronList), vertices(particleList)
{}

void Sim::init()
{
	// Precompute rest deformation (Dm), volume, inverse Dm for each tetrahedron
	for(int i=0; i<tetras->numTetra; i++)
	{
		tetras->computeRestDeformation( i, vertices );
		tetras->computeInvRestDeformation( i );
		tetras->computeUndeformedVolume( i );
	}
}

void Sim::clean()
{
	//Set forces for all vertices/particles to zero
	for(int i=0; i<vertices->numParticles; i++)
	{
		vertices->force[i] = Eigen::Matrix<T, 3, 1>::Zero();
	}
}

void Sim::eulerIntegration()
{
// TODO
}

void Sim::computeElasticForces( int tetraIndex )
{
	Eigen::Matrix<T,3,3> newDeformation = tetras->computeNewDeformation( tetraIndex, vertices ); // Compute Ds, the new deformation
	Eigen::Matrix<T,3,3> F 				= tetras->computeF( tetraIndex, newDeformation ); // Compute F = Ds(Dm_inv)
	Eigen::Matrix<T,3,3> P 				= tetras->computeP( tetraIndex, F ); // Compute Piola (P)
	Eigen::Matrix<T,3,3> H 				= tetras->computeH( tetraIndex, P, newDeformation ); // Compute Energy (H)

	tetras->addForces( tetraIndex, vertices, H );// Add energy to forces (f += h)
}

void Sim::update()
{
	clean(); //clears forces
	checkCollisions(); //Apply Forces to particles that occur through collision

	// Loop through tetras
	for(int i=0; i<tetras->numTetra; i++)
	{
		computeElasticForces(i); //computes and adds elastic forces to each particle
	}

	for(int i=0; i<tetras->numTetra; i++)
	{
		eulerIntegration();
	}
}

void Sim::checkCollisions()
{
	// TODO
	// First do brute force SDF for primitives
	// Later do interobject collisions with a grid+bounding box or BVH

	for(int i=0; i<vertices->numParticles; i++)
	{
		//Apply force to particle once not per tetrahedron
		//use hash table for this
	}
}
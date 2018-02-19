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

void Sim::computeElasticForce( int tetraIndex )
{
	Eigen::Matrix<T, 3, 3>::Zero(); //Ds

	// Compute Ds, the new deformation
	// Eigen::Matrix<T,3,3> newDeformation = tetras->computeNewDeformation( tetraIndex, vertices );

	// Eigen::Matrix<T,3,3> F = tetras->computeF( tetraIndex ); // Compute F = Ds(Dm inv)
	// Eigen::Matrix<T,3,3> P = tetras->computeP( tetraIndex, vertices ); // Compute Piola (P)
	// Eigen::Matrix<T,3,3> H = tetras->computeH( tetraIndex, vertices ); // Compute Energy (H)

	// Add energy to forces (f += h)
	// f4 += -(h1 + h2 + h3)
}

void Sim::update()
{
	clean(); //clears forces
	checkCollisions(); //Apply Forces to particles that occur through collision

	// Loop through tetras
	for(int i=0; i<tetras->numTetra; i++)
	{
		computeElasticForce(i);
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

// Reads all the required files and stores stuff in respective arrays
//void tetRead()
//{
//    // read node file and store list of points (id, pos, attributes)
//    // read .ele file and store tetrahedra (id, nodes)
//    // can read more files for faces, edges, etc.
//}
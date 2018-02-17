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

void Sim::eulerIntegration()
{
    // TODO
}

void Sim::update()
{
    // TODO
    // ***** BETTER DOUBLE CHECK WITH LADISLAV'S VIDEO BOI *****

    // Loop through all the particles
    //     forces = 0
    //     checkCollision()
    //          Apply force to particle once not tetrahedron 

    // Loop through tetras
    //      Compute new deformation (Ds)
    //      Compute force = Ds(Dm inv)
    //      Compute Piola (P) -- Need a table of values
    //      Compute Energy (H)
    //      Add energy to forces (f += h)
    //      f4 += -(h1 + h2 + h3)
    //      eulerIntegration()
}

void Sim::checkCollisions()
{
    // TODO
    // First do brute force SDF for primitives
    // Later do interobject collisions with a grid+bounding box or BVH
}

// Reads all the required files and stores stuff in respective arrays
//void tetRead()
//{
//    // read node file and store list of points (id, pos, attributes)
//    // read .ele file and store tetrahedra (id, nodes)
//    // can read more files for faces, edges, etc.
//}
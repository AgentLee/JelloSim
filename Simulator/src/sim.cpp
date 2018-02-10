/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

void init()
{
    // TODO
    // Precompute rest deformation (Dm), volume, inverse Dm
}

void eulerIntegration()
{
    // TODO
}

void update()
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

void checkCollisions()
{
    // TODO
    // First do brute force SDF for primitives
    // Later do interobject collisions with a grid+bounding box or BVH
}

void sim()
{
    std::vector<Tetrahedron> tetras;    
    std::shared_ptr<Particles> particles;
}
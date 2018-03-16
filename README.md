Jello Simulator using FEM
=========================

## Overview
This project implements a Jello Simulator using the Finite Element Method (FEM). The Finite Element method approximates the values of the unknowns (in this case, forces) at discrete points over the domain of the simulation. It does this by breaking down the problem into a number of really tiny simple problems, called finite elements. We treat these finite elements as their own systems that we solve individually to obtain a solution for the overall system. Because the finite element method uses approximations of the system at the discrete points its accuracy increase with smaller and smaller elements.

The simulator can actually represent a whole host of elastic materials by specifing their Densities, Young's Modulus and Poisson's Ratio; We have however set it up with values that represent Jello.

Our simulator is set up to create files (BGEO files) that can be rendered out in Houdini.

Refer to the file labeled "USAGE_INSTRUCTIONS" for how to use the Simulator.

#### Contributors
- [Aman Sachan](https://github.com/aman-sachan-asach)
- [Jonathan Lee](https://github.com/agentlee)
- [Rishabh Shah](https://github.com/rms13)
- [Grover Chen](https://github.com/Groverc85)

#### Skip Forward to:
1. [Features](#Features)
2. [Implementation Overview](#Implementation)
3. [Algorithm Layout](#Algorithm)
4. [Architecture](#Architecture)
5. [Future Work](#Future)
6. [Limitations](#Limitations)
7. [Resources](#Resources)

## Features <a name="Features"></a>
 - Finite Element Method with Fixed Corotated Elastic Model
 - Collisions with solid objects
 - Fixed point constraints
 - Scene creation and Rendering
 - Mass Distribution across discretised mesh
 - Explicit Forward Euler integration scheme
 - Data Driven Architecture that's easy to understand and extend
 - BGEO writing facilities
 - Obj to Poly file conversion utilities to make tetgen viable for meshes
 - Inter-Mesh Collisions

## Implementation Overview <a name="Implementation"></a>

The simulation works by applying forces to the individual vertices of that make up the mesh. These vertices exist through out the entire volume of the mesh and not just on it's surface. The forces are then used to update the velocities and then positions of the vertices in the Euler Integration Step. Because the vertices have changed positions the triangles that represent the surface of the mesh have also changed in size, shape, and position. This is what is seen as the elastic deformation that occurs on the model.

The forces that are applied to the mesh are both internal and external. The _external forces_ acting on a mesh are things like gravity or wind forces. These forces are applied in the same magnitude to all of the vertices that are exposed to those forces.
_Internal forces_ are the more interesting forces, as these are responsible for the elastic behaviour of the mesh. Basically, we assume that whatever shape the mesh starts off in is its default shape free of any deformation. And when the vertices change positions due to collisions, the changes in position are not uniform across the mesh. This results in deformation, a change of form of the mesh. We convert this deformation into forces that want to change the mesh back to an undeformed state.

So, far we have talked about the Mesh as a single system, but because we are using the Finite Element Method, we actually represent it as a bunch of tiny tetrahedrons, that are the finite elements. These tetrahedrons are the ones that have an initial undeformed state and undergo deformations as described above. The cummulation of all the tetrahedrons solves the bigger system that is the mesh.
The finer the tetrahedrons that make up the mesh the slower but also more accurate the resulting simulation will be.

### Algorithm Layout <a name="Algorithm"></a>

The algorithm can be broken down into a _Preprocess Step_ and an _Update Loop_.

#### Preprocess Step/Initialization

```
For every Mesh in the Simulation
{
   For every Tetrahedron that makes up that Mesh
   {
      compute the Rest Deformation for that tetrahedron
      compute the Inverse of the Rest Deformation for that tetrahedron
      compute the Undeformed Volume for that tetrahedron
      compute the product of (Undeformed Volume * Inverse Rest Deformation) for that tetrahedron --> Optimization Purposes
      
      compute the Distribution of Mass for that tetrahedron
   }
}
```

#### Simulation Update Loop

```
Update
{
   Clear Forces on all Vertices
   Recompute Mesh Attributes like triangle normals and bounding boxes
   Compute and Add External Forces for all Vertices
   Compute and Add Internal Forces for all Vertices
   Euler Integration For every Mesh in the Sim
   Collision Handling For every Mesh in the Sim
}
```

#### Elastic Force Calculation Loop

```
For every Mesh in the Simulation
{
   For every Tetrahedron that makes up that Mesh
   {
      Compute the New Deformation of that tetrahedron
      Compute the F Matrix for that tetrahedron
      Compute the P Matrix for that tetrahedron
      Compute the H Matrix for that tetrahedron
      
      The individual columns of the H Matrix represent the forces dues to deformation on the first 3 vertices of the tetrahedron.
      The Force on the fourth vertex = -(sum of forces on the other 3 vertices)
      
      Add forces to the individual vertices of that tetrahedron
   }
}
```

F can be thought of as the Transformation Matrix that converts the Undeformed tetrahedron into the Deformed Tetrahedron.
P is the Piola Kirchoff Stress Tensor that is a force computation model.
H is a matrix of Forces on the individual vertices.

### Collisions with rigid objects:
### Fixed Point Constraints:


## Architecture <a name="Architecture"></a>

## Future Work <a name="Future"></a>
 - Parallelization through CUDA so that you can run it on the GPU

## Limitations <a name="Limitations"></a>
 - Not energy conserving -- the mesh does not keep bouncing forever
 - Physics gets fudged a little during collisions (both with other ‘jello’ meshes and with sdf’s) - we don’t calculate momentum transfer but rely on deformation, velocity, and position changes to fix most things.

## Resources <a name="Resources"></a>
- http://run.usc.edu/femdefo/sifakis-courseNotes-TheoryAndDiscretization.pdf
- http://run.usc.edu/femdefo/barbic-courseNotes-modelReduction.pdf
- http://www.femdefo.org/
- http://www.math.ucla.edu/~jteran/papers/ITF04.pdf
- https://www.youtube.com/watch?v=BH1OrCtaPjo
- https://www.youtube.com/watch?v=SACGMSZx4FY&t=1s

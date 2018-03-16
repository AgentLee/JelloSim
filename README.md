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
3. [Algorithm Pseudo Code](#Algorithm)
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

### Algorithm Pseudo Code <a name="Algorithm"></a>
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

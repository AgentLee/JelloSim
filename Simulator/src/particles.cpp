#include "particles.h"

using T = double;

Particles::Particles(int n, T initialMass): numParticles(n)
{
    for(int i=0; i<numParticles; i++)
    {
        mass.push_back(initialMass);
        
        Eigen::Matrix<T,1,3> vel_particle = Eigen::Matrix<T,1,3>::Zero();
        vel.push_back(vel_particle);

        Eigen::Matrix<T,1,3> pos_particle = Eigen::Matrix<T,1,3>::Zero();
        pos.push_back(pos_particle);

        Eigen::Matrix<T,1,3> force_particle = Eigen::Matrix<T,1,3>::Zero();
        force.push_back(force_particle);
    }
}

/* 
 * Initialize all particles with random values
 */
void Particles::initializeParticles_random()
{
    for(int i=0; i<numParticles; i++)
    {
        mass[i] = std::rand() % 100;
        vel[i] = Eigen::Matrix<T,1,3>::Random(1,3);
        pos[i] = Eigen::Matrix<T,1,3>::Random(1,3);
    }
}

/* 
 * Using sympledic Euler to update speed and velocity based on 
 * the force, mass and previous state
 */
void Particles::updateParticlePositions(T dt)
{
    for(int i=0; i<numParticles; i++)
    {
        pos[i] += vel[i] * dt;
    }
}

void Particles::updateParticleVelocity(T dt)
{
    for(int i=0; i<numParticles; i++)
    {
        vel[i] += force[i] / mass[i] * dt;
    }
}

/*
 *  .node FILE FORMAT:
 *      http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html
 *
 *  First line:
 *              <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
 *
 *  Remaining lines list # of points:
 *              <point #> <x> <y> <z> [attributes] [boundary marker]
 *
 */

// helper function
void Particles::tetgen_readLine(std::ifstream &fin, int numDims, int numAtt)
{
    float f;
    fin >> f; // first one is id..
    for(int i = 0; i < numDims; ++i)
    {
        fin >> pos[f-1](0,i);
    }
    // std::cout << pos[f-1](0, 0) << " " << pos[f-1](1, 0) << " " << pos[f-1](2, 0) << "\n";

    for(int i = 0; i < numAtt; ++i)
    {
        float someAttribute;
        fin >> someAttribute;
    }
}

void Particles::tetgen_readNode(const std::string &inputFileName)
{
    // TODO
    std::ifstream fin(inputFileName);

    if(fin.is_open())
    {
        int numPoints;
        int numDims;
        int numAtts;
        int boundaryMarker;

        fin >> numPoints >> numDims >> numAtts >> boundaryMarker;

        //Resizing all vectors in Particles
        numParticles = numPoints;
        mass.resize(numPoints, 1.0f); //giving mass of 1 to all particles
        vel.resize(numPoints, Eigen::Matrix<T,1,3>::Zero());
        pos.resize(numPoints);
        force.resize(numPoints);

        for(int i = 0; i < numPoints; ++i)
        {
            tetgen_readLine(fin, numDims, numAtts);
        }

        fin.close();
    }
}

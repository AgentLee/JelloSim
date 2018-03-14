#include "particles.h"
#include <Partio.h>

Particles::Particles(int n, T initialMass): numParticles(n)
{
    for(int i=0; i<numParticles; i++)
    {
        mass.push_back(initialMass);
        
        Eigen::Matrix<T, 3, 1> vel_particle = Eigen::Matrix<T, 3, 1>::Zero();
        vel.push_back(vel_particle);

        Eigen::Matrix<T, 3, 1> pos_particle = Eigen::Matrix<T, 3, 1>::Zero();
        pos.push_back(pos_particle);

        Eigen::Matrix<T, 3, 1> force_particle = Eigen::Matrix<T, 3, 1>::Zero();
        force.push_back(force_particle);
    }
}

/* 
 * Using sympledic Euler to update speed and velocity based on 
 * the force, mass and previous state
 */
void Particles::updateAllParticlePositions(T dt)
{
    for(int i=0; i<numParticles; i++)
    {
        pos[i] += vel[i] * dt;
    }
}

void Particles::updateAllParticleVelocities(T dt, Bounds& FixedRegion)
{
    for(int i=0; i<numParticles; i++)
    {
        Point3f position;
        position << pos[i][0], pos[i][1], pos[i][2];
        if(FixedRegion.Contains(position))
        {
            vel[i] *= 0.0;
        }
        else
        {
            vel[i] += force[i] / mass[i] * dt;
        }
    }
}

void Particles::updateParticlePosition(T dt, uint index)
{
    pos[index] += vel[index] * dt;
}

void Particles::updateParticleVelocity(T dt, uint index)
{
    vel[index] += force[index] / mass[index] * dt;
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
        fin >> pos[f-1](i);
    }

    for(int i = 0; i < numAtt; ++i)
    {
        float someAttribute;
        fin >> someAttribute;
    }
}

void Particles::tetgen_readNode(const std::string &inputFileName)
{
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
        mass.resize(numPoints, 0.0f); //giving mass of 1 to all particles
        vel.resize(numPoints, Eigen::Matrix<T, 3, 1>::Zero());
        pos.resize(numPoints);
        force.resize(numPoints);

        for(int i = 0; i < numPoints; ++i)
        {
            tetgen_readLine(fin, numDims, numAtts);
        }

        fin.close();
    }
}
#include "particles.h"
#include <Partio.h>

using namespace std;
using namespace glm;


void SolveGD(mat3 A, vec3 & x, vec3 b) {
    
    // NOTE, operators are overloaded by LargeVM
    
    vec3 r = b - (A * x);
    vec3 q;
    
    float alpha;
   
    float del = dot(r, r);
    float del0 = del;
    
    float eTol = std::pow(1*10, -10); // Error tolerance
    int iterMax = 10; // maximum iterations
    int i = 0;
    
    while ((del > (std::pow(eTol, 2) * del0)) && (i < iterMax)) {
        
        q = A * r;
        
        alpha = del / dot(r, q);
        
        x = x + alpha * r;
        r = r - alpha * q;
        
        del = dot(r, r);
        
        i++;
    }
}


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

void Particles::updateAllParticleVelocities(T dt)
{
    for(int i=0; i<numParticles; i++)
    {
        vel[i] += force[i] / mass[i] * dt;
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


/**
* Simulate.
* Note external forces must have been computed prior to invocation.
*
* The dynamic equation has the following form (where u=x-x_0, and x' is derivative wrt. time)
*
*   M x'' + Cx' + K (x-x_0) = f_ext
*
* This can be transformed to a system of 2x3n equations of first order derivative:
*
*     x' = v
*   M v' = - C v - K (x-x_0) + f_ext
*
* The semi-implicit Euler scheme approximates the above with:
*
*     x^(i+1) = x^i + \delta t * v^(i+1)
*   M v^(i+1) = M v^i + \delta t ( - C v^(i+1)  - K ( x^i - x_0 ) + f^i_ext
*
* This is solved using implicit integration.
*/


void Particles::computeForce(int index) {
    glm::vec3 deltaP = glm::vec3((pos[index+1] - pos[index])(0), (pos[index+1] - pos[index])(1), (pos[index+1] - pos[index])(2));
    glm::vec3 deltaP2 = deltaP * deltaP;
    glm::vec3 pdC_pdX = deltaP / (float)1.0;

    double C = 1.0;
    glm::vec3 V = glm::vec3((float) vel[index](0), (float) vel[index](1), (float) vel[index](2));
    double C_DOT = glm::dot(V, -pdC_pdX) + glm::dot(V, pdC_pdX);

    glm::mat3 pdF_pdX, pdF_pdV;
    glm::mat3 pd2C_pdX2[2][2];
    
    pd2C_pdX2[0][0][0][0] = (-C * deltaP2[0]) + C;
    pd2C_pdX2[0][0][1][1] = (-C * deltaP2[1]) + C;
    pd2C_pdX2[0][0][2][2] = (-C * deltaP2[2]) + C;
    
    pd2C_pdX2[0][1][0][0] = (C * deltaP2[0]) - C;
    pd2C_pdX2[0][1][1][1] = (C * deltaP2[1]) - C;
    pd2C_pdX2[0][1][2][2] = (C * deltaP2[2]) - C;
    
    pd2C_pdX2[1][0] = pd2C_pdX2[0][1];
    pd2C_pdX2[1][1] = pd2C_pdX2[0][0];
    
    glm::mat3 pdX1 = outerProduct(pdC_pdX, pdC_pdX);
    glm::mat3 pdX2 = outerProduct(pdC_pdX, -pdC_pdX);
    glm::mat3 pdX3 = outerProduct(-pdC_pdX, -pdC_pdX);
    
    pdF_pdX += - ((float)Ks * (pdX1 + (pd2C_pdX2[0][0] * (float)C))) - ((float)Kd * (pd2C_pdX2[0][0] * (float)C_DOT));
    pdF_pdX += - ((float)Ks * (pdX2 + (pd2C_pdX2[0][1] * (float)C))) - ((float)Kd * (pd2C_pdX2[0][1] * (float)C_DOT));
    pdF_pdX += - ((float)Ks * (pdX3 + (pd2C_pdX2[1][1] * (float)C))) - ((float)Kd * (pd2C_pdX2[1][1] * (float)C_DOT));
    
    pdF_pdV += - (float)Kd * pdX1;
    pdF_pdV += - (float)Kd * pdX2;
    pdF_pdV += - (float)Kd * pdX3;
    
    dForceDX.push_back(pdF_pdX);
    dForceDV.push_back(pdF_pdV);
}

void Particles::updateAllParticlePositionsBackwardEuler(T dt) {
    for(int i=0; i<numParticles; i++)
    {
        computeForce(i);

        glm::vec3 Vnew;

        glm::mat3 A = mMat - ((float)(pow(dt, 2))*kMat);
        glm::vec3 F = glm::vec3((float) force[i](0), (float) force[i](1), (float) force[i](2));
        glm::vec3 b = (mMat* V) + ((float)dt * F);

        SolveGD(A, Vnew, b);

        pos[i](0) += dt * Vnew[i];
    }
}

void Particles::updateAllParticleVelocitiesBackwardEuler(T dt) {
    for(int i=0; i<numParticles; i++)
    {
        computeForce(i);

        glm::vec3 Vnew;

        glm::mat3 A = mMat - ((float)(pow(dt, 2))*kMat);
        glm::vec3 F = glm::vec3((float) force[i](0), (float) force[i](1), (float) force[i](2));
        glm::vec3 b = (mMat* V) + ((float)dt * F);

        SolveGD(A, Vnew, b);

        vel[index] = Vnew[index];
    }
}


void Particles::updateParticlePositionBackwardEuler(T dt, uint index) {
    computeForce(i);

    glm::vec3 Vnew;

    glm::mat3 A = mMat - ((float)(pow(dt, 2))*kMat);
    glm::vec3 F = glm::vec3((float) force[i](0), (float) force[i](1), (float) force[i](2));
    glm::vec3 b = (mMat* V) + ((float)dt * F);

    SolveGD(A, Vnew, b);

    pos[index] += dt * Vnew[index];
}

void Particles::updateParticleVelocityBackwardEuler(T dt, uint index) {
    computeForce(i);

    glm::vec3 Vnew;

    glm::mat3 A = mMat - ((float)(pow(dt, 2))*kMat);
    glm::vec3 F = glm::vec3((float) force[i](0), (float) force[i](1), (float) force[i](2));
    glm::vec3 b = (mMat* V) + ((float)dt * F);

    SolveGD(A, Vnew, b);

    vel[index] = Vnew[index];
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

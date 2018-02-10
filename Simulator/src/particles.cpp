#include "particles.h"

using T = double;

Particles::Particles(int n, T initialMass): numParticles(n)
{
	for(int i=0; i<numParticles; i++)
	{
		mass.push_back(initialMass);
		
		Eigen::Matrix<T,3,1> vel_particle = Eigen::Matrix<T,3,1>::Zero();
		vel.push_back(vel_particle);

		Eigen::Matrix<T,3,1> pos_particle = Eigen::Matrix<T,3,1>::Zero();
		pos.push_back(pos_particle);
	}
}

void Particles::initializeParticles_random()
{
	for(int i=0; i<numParticles; i++)
	{
		mass[i] = std::rand() % 100;
		vel[i] = Eigen::Matrix<T,3,1>::Random(3,1);
		pos[i] = Eigen::Matrix<T,3,1>::Random(3,1);
	}
}

void Particles::updateParticlePositions(T dt)
{
	for(int i=0; i<numParticles; i++)
	{
		pos[i](0,0) += vel[i](0,0)*dt;
		pos[i](1,0) += vel[i](1,0)*dt;
		pos[i](2,0) += vel[i](2,0)*dt;
	}
}

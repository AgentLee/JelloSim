/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SET_POSITIONS 1
#define SET_VELOCITIES 0
#define PAPER 0
#define INTER_OBJECT_COLLISIONS 1

Sim::Sim( std::vector<std::shared_ptr<Mesh>>& MeshList ) : MeshList(MeshList)
{}

void Sim::init()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

		// Precompute rest deformation (Dm), volume, inverse Dm, and volume*inverseDmTranspose for each tetrahedron
		tetras->restDeformation.resize(tetras->numTetra);
		tetras->restInverseDeformation.resize(tetras->numTetra);
		tetras->undeformedVolume.resize(tetras->numTetra);
		tetras->undefVol_into_restInvDefTranspose.resize(tetras->numTetra);

		for(int i=0; i<tetras->numTetra; i++)
		{
			tetras->computeRestDeformation( i, vertices );
			tetras->computeInvRestDeformation( i );
			tetras->computeUndeformedVolume( i );
			tetras->computeUndefVol_into_restInvDefTranspose( i );

			tetras->addMass( i, vertices );
		}

        std::cout << MeshList[j]->AABB.min[0] << " " << MeshList[j]->AABB.min[1] << " " << MeshList[j]->AABB.min[2] << std::endl;
        std::cout << MeshList[j]->AABB.max[0] << " " << MeshList[j]->AABB.max[1] << " " << MeshList[j]->AABB.max[2] << std::endl;
	}

}

void Sim::clean()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
		//Set forces for all vertices/particles to zero
		for( int i=0; i < vertices->numParticles; i++ )
		{
			vertices->force[i] = Eigen::Matrix<T, 3, 1>::Zero();
		}
	}
}

void Sim::eulerIntegration(float dt)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
        SDF_Collisions(dt, i);
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		vertices->updateAllParticleVelocities(dt);
		vertices->updateAllParticlePositions(dt);
	}
}

void Sim::eulerIntegrationWithCollisionTesting(float dt)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		SDF_Collisions(dt, i);

		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		vertices->updateAllParticleVelocities(dt);
		vertices->updateAllParticlePositions(dt);
	}

	for(uint j=0; j<MeshList.size(); j++)
	{
		//Check if this mesh's AABB is intersecting with any other mesh's AABB
		for(uint i=0; i<MeshList.size() && i!=j; i++)
		{
			bool intersects = Intersect_AABB_with_AABB( MeshList[j]->AABB, MeshList[i]->AABB );
			
			if(intersects)
			{
				std::shared_ptr<Particles> vertices = MeshList[j]->vertices; //Vertices of Current Mesh
				// Then do a brute force check, i.e loop through all vertice of one mesh 
				// OR use a grid structure or cull the triangles somehow
				
				for(int k = 0; k< vertices->numParticles; ++k)
				{
					bool collided = false;
					Mesh_Collisions(dt, j, i, vertices, k, collided);

					if(!collided)
					{
						vertices->updateParticleVelocity(dt, k);
						vertices->updateParticlePosition(dt, k);
					}

                    Point3f p;
                    p << vertices->pos[k][0], vertices->pos[k][1], vertices->pos[k][2];
				}
			}
		}
	}
}

void Sim::addExternalForces()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;
		
		for(int i=0; i<vertices->numParticles; i++)
		{
			vertices->force[i](1) -= 9.81 * vertices->mass[i]; // gravity
		}
	}
}

void Sim::computeElasticForces( int frame )
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
		std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

		// Loop through tetras
		for(int tetraIndex=0; tetraIndex < tetras->numTetra; tetraIndex++)
		{
			Eigen::Matrix<T,3,3> newDeformation = tetras->computeNewDeformation( tetraIndex, vertices ); // Compute Ds, the new deformation
			Eigen::Matrix<T,3,3> F 				= tetras->computeF( tetraIndex, newDeformation ); // Compute F = Ds(Dm_inv)
			Eigen::Matrix<T,3,3> P 				= tetras->computeP( tetraIndex, F, frame ); // Compute Piola (P)
			Eigen::Matrix<T,3,3> H 				= tetras->computeH( tetraIndex, P ); // Compute Energy (H)

			tetras->addForces( tetraIndex, vertices, H );// Add energy to forces (f += h)
		}
	}
}

void Sim::update(float dt, int frame)
{
	clean(); //clears forces
    reComputeMeshAttributes(); //compute triangle normals for all meshes

    addExternalForces();
	computeElasticForces(frame); //computes and adds elastic forces to each particle
#if INTER_OBJECT_COLLISIONS
    eulerIntegrationWithCollisionTesting(dt);
#else
    eulerIntegration(dt);
#endif
}

void Sim::reComputeMeshAttributes()
{
	for(uint j=0; j<MeshList.size(); j++)
	{
		MeshList[j]->triangles->computeNormals(MeshList[j]->vertices);
        MeshList[j]->calcBounds();
	}
}

// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
bool Sim::LineTriangleIntersection(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, Intersection *isect)
{
	// create ray and call triangle intersection for all triangles..
	// store triangle reference
	Ray r;
	r.origin = origPos;
	r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);
    float length = r.direction.norm();
    r.direction.normalize();

	// assuming only 2 meshes for now..
	int tris = MeshList[1]->triangles->triFaceList.size();
	isect->t = -1;
    isect->hit = false;
	isect->triangleIndex = -1;
	for(int i=0; i < tris; i++)
	{
		float tTemp = isect->t;
        Eigen::Matrix<T, 3, 1> baryCoords;
		bool intersect = MeshList[1]->triangles->intersect(r, i, &tTemp, MeshList[1]->vertices, &baryCoords);
		if(intersect && isect->t > tTemp && tTemp <= length)
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp;
            isect->BarycentricWeights = baryCoords;
            isect->normal = MeshList[1]->triangles->triNormalList[i];
            isect->point = r.origin + tTemp * r.direction;
		}
	}

	return isect->hit;
}

void Sim::resolveCollisions( std::shared_ptr<Triangles>& triangles, std::shared_ptr<Particles>& vertices, 
	Intersection& isect, Eigen::Matrix<T, 3, 1>& displacement,
	Eigen::Matrix<T, 3, 1>& particlePos, Eigen::Matrix<T, 3, 1>& particleVel )
{
	// NOTE: triangles and vertices corresspond to the Triangles and Vertices of the other Mesh

	// set position or velocity or both OR use momentum OR use paper's implementation with normal reaction forces and friction
	Eigen::Matrix<uint, 3, 1> verticesOfTriangle;
    verticesOfTriangle = triangles->triFaceList[isect.triangleIndex];

	if(SET_POSITIONS)
	{
		//---------------- Setting Positions ---------------------------------
		particlePos = isect.point.cast<T>() - 0.25f*displacement; // 1/4th to moving vertex
		displacement = 0.75f*displacement;				// 3/4th to moving triangle

		// move vertices of triangle according to barycentric weights
		vertices->pos[verticesOfTriangle[0]] += isect.BarycentricWeights[0] * displacement;
		vertices->pos[verticesOfTriangle[1]] += isect.BarycentricWeights[1] * displacement;
		vertices->pos[verticesOfTriangle[2]] += isect.BarycentricWeights[2] * displacement;
	}
	if(SET_VELOCITIES)
	{
		particleVel = Eigen::Matrix<T, 3, 1>::Zero();

		vertices->vel[verticesOfTriangle[0]] = Eigen::Matrix<T, 3, 1>::Zero();
		vertices->vel[verticesOfTriangle[1]] = Eigen::Matrix<T, 3, 1>::Zero();
		vertices->vel[verticesOfTriangle[2]] = Eigen::Matrix<T, 3, 1>::Zero();
	}
	if(PAPER)
	{
	}
}

void Sim::SDF_Collisions(float dt, uint j)
{
	// Check if the Mesh hit the ground or any other solid piece of geometry that is 
	// essentially an infinite Mass Rigid Body that doesnt move
	std::shared_ptr<Tetrahedrons> tetras = MeshList[j]->tetras;
	std::shared_ptr<Particles> vertices = MeshList[j]->vertices;

	for(int i = 0; i< vertices->numParticles; ++i)
	{
		Eigen::Matrix<T, 3, 1> p(vertices->pos[i]);

		// Transform the vertex to the plane's local space
		// Assume the plane is at the origin
		Eigen::Matrix<T, 4, 1> n = Eigen::Matrix<T, 4, 1>(0, 1, 0, 0);
		float sdf = SDF::sdPlane(p, n);

		// Check if particle went through the surface
		if(sdf < 0) 
		{
			if( vertices->vel[i].dot(Eigen::Matrix<T, 3, 1>(0, 1, 0)) < 0 )
			{
				vertices->vel[i] = Eigen::Matrix<T, 3, 1>::Zero();
			}
		}
	}
}

void Sim::Mesh_Collisions(float dt, uint i, uint j, std::shared_ptr<Particles>& vertices, int particleIndex, bool& collided)
{
	std::shared_ptr<Triangles> triangles = MeshList[j]->triangles;

	// for every vertex create a ray from the current position to its projected position in the next frame
	// See if that ray intersects any triangle that Belongs to the other mesh.
	Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[particleIndex] + vertices->vel[particleIndex] * dt;
	Intersection isect;
	LineTriangleIntersection(vertices->pos[particleIndex], projectedPos, &isect);

	if(isect.hit)
	{
		//resolve collisions between vertex and a triangle
		Eigen::Matrix<T, 3, 1> displscement = projectedPos - vertices->pos[j];

		resolveCollisions( MeshList[i]->triangles, MeshList[i]->vertices, isect, displscement, 
						vertices->pos[particleIndex], vertices->vel[particleIndex] );
		collided = true;
	}
}

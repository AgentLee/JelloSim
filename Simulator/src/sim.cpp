/*
 * This file contains the main simulation loop for FEM.
 */

#include "sim.h"

#define SET_POSITIONS 1
#define SET_VELOCITIES 1
#define PAPER 0
#define INTER_OBJECT_COLLISIONS 1

Sim::Sim( std::vector<std::shared_ptr<Mesh>>& MeshList, Bounds& FixedRegion ) : MeshList(MeshList), FixedRegion(FixedRegion)
{}

void Sim::init()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->initMeshForSim();
	}
}

void Sim::clean()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->clearForces();
	}
}

void Sim::reComputeMeshAttributes()
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		MeshList[i]->triangles->computeNormals(MeshList[i]->vertices);
		MeshList[i]->calcBounds();
	}
}

void Sim::computeAllForces( int frame )
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		//////////------------------ EXTERNAL FORCES ---------------/////////////
		//Not abstracting to mesh because we could add forces that only operate on 
		//some of the vertices and things like that which is annoying to genarlize for
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		
		for(int j=0; j<vertices->numParticles; j++)
		{
			vertices->force[j](1) -= 9.81 * vertices->mass[j]; // gravity
		}

		//////////------------------ ELASTIC FORCES ---------------/////////////
		MeshList[i]->computeElasticForcesOnMesh(frame);
	}
}

void Sim::update(float dt, int frame)
{
	clean(); //clears forces
    reComputeMeshAttributes(); //compute triangle normals for all meshes
	computeAllForces(frame); //computes and adds elastic forces to each particle

#if INTER_OBJECT_COLLISIONS
    eulerIntegration(dt);
    Collisions(dt);
#else
    eulerIntegrationWithSDF_Collisions(dt);
#endif
}

void Sim::eulerIntegration(float dt)
{
    for(uint i=0; i<MeshList.size(); i++)
    {
        MeshList[i]->setPrevVerts();
        std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
        vertices->updateAllParticleVelocities(dt, FixedRegion);
        vertices->updateAllParticlePositions(dt);
    }
}

void Sim::eulerIntegrationWithSDF_Collisions(float dt)
{
	for(uint i=0; i<MeshList.size(); i++)
	{
        SDF_Collisions(dt, i);
		std::shared_ptr<Particles> vertices = MeshList[i]->vertices;
		vertices->updateAllParticleVelocities(dt, FixedRegion);
		vertices->updateAllParticlePositions(dt);
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

void Sim::Collisions(float dt)
{
    for(uint i=0; i<MeshList.size(); i++)
    {
        SDF_Collisions(dt, i);
    }

    std::vector<std::vector<Eigen::Matrix<T, 3, 1>>> newPositions;
    std::vector<std::vector<Eigen::Matrix<T, 3, 1>>> newVelocities;
    for(uint i = 0; i < MeshList.size(); i++)
    {
        std::vector<Eigen::Matrix<T, 3, 1>> pos;
        std::vector<Eigen::Matrix<T, 3, 1>> vel;
        for(int j = 0; j < MeshList[i]->vertices->numParticles; j++)
        {
            pos.push_back(MeshList[i]->vertices->pos[j]);
            vel.push_back(MeshList[i]->vertices->vel[j]);
        }
        newPositions.push_back(pos);
        newVelocities.push_back(vel);
    }


    for(uint i=0; i<MeshList.size(); i++)
    {
        for (uint j =i+1; j < MeshList.size(); j++)
        {
            std::shared_ptr<Mesh> fallingMesh = MeshList[i];
            std::shared_ptr<Mesh> standingMesh = MeshList[j];

            bool possibleIntersectsion = Intersect_AABB_with_AABB(fallingMesh->AABB, standingMesh->AABB);

            if (possibleIntersectsion) 
            {
                for (int k = 0; k < fallingMesh->vertices->numParticles; ++k) 
                {
                    Intersection isect;
                    IntersectionTesting( fallingMesh, standingMesh, k, &isect);

                    if (isect.hit) 
                    {
                        std::shared_ptr<Triangles> triangles = standingMesh->triangles;

                        //resolve collisions between vertex and a triangle
                        if (SET_POSITIONS)
                        {
                            newPositions[i][k] = fallingMesh->prevVerts->pos[k];
                            // newPositions[i][k] = isect.point - 2.0 * isect.penetrationDistance;

                            // move vertices of triangle according to barycentric weights
                            newPositions[j][isect.vertsOfTriangle[0]] += isect.penetrationDistance;
                            newPositions[j][isect.vertsOfTriangle[1]] += isect.penetrationDistance;
                            newPositions[j][isect.vertsOfTriangle[2]] += isect.penetrationDistance;
                            newPositions[j][isect.vertsOfTriangle[0]] = standingMesh->prevVerts->pos[isect.vertsOfTriangle[0]];
                            newPositions[j][isect.vertsOfTriangle[1]] = standingMesh->prevVerts->pos[isect.vertsOfTriangle[1]];
                            newPositions[j][isect.vertsOfTriangle[2]] = standingMesh->prevVerts->pos[isect.vertsOfTriangle[2]];
                        }
                        if (SET_VELOCITIES)
                        {
                            // Eigen::Matrix<T, 3, 1> normal;
                            // normal = triangles->triNormalList[isect.triangleIndex];

                            // Eigen::Matrix<T, 3, 1> d;
                            // d = fallingMesh->vertices->vel[k];

                            // if(normal.dot(d) < 0.0f)
                            // {
                            //     normal = -normal;
                            // }

                            // d.normalize();
                            // Eigen::Matrix<T, 3, 1> reflected;
                            // reflected = d - 2.0f*(d.dot(normal))*normal;
                            // newVelocities[i][k] = reflected;//Eigen::Matrix<T, 3, 1>::Zero();
                            newVelocities[i][k] = Eigen::Matrix<T, 3, 1>::Zero();

                            newVelocities[j][isect.vertsOfTriangle[0]] += fallingMesh->vertices->vel[k];
                            newVelocities[j][isect.vertsOfTriangle[1]] += fallingMesh->vertices->vel[k];
                            newVelocities[j][isect.vertsOfTriangle[2]] += fallingMesh->vertices->vel[k];
                        }

                        fallingMesh->vertices->color[k] = isect.color;
                        standingMesh->vertices->color[isect.vertsOfTriangle[0]] = isect.color;
                        standingMesh->vertices->color[isect.vertsOfTriangle[1]] = isect.color;
                        standingMesh->vertices->color[isect.vertsOfTriangle[2]] = isect.color;
                    }
                }
            }
        }
    }

    for(uint i = 0; i < MeshList.size(); i++)
    {
        for(int j = 0; j < MeshList[i]->vertices->numParticles; j++)
        {
            MeshList[i]->vertices->pos[j] = newPositions[i][j];
            MeshList[i]->vertices->vel[j] = newVelocities[i][j];
        }
    }
}

// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
bool Sim::IntersectionTesting(std::shared_ptr<Mesh> meshA, std::shared_ptr<Mesh> meshB, int vertexID, Intersection *isect)
{
    Eigen::Matrix<T, 3, 1>& vCurr = meshA->vertices->pos[vertexID];
    Eigen::Matrix<T, 3, 1>& vPrev = meshA->prevVerts->pos[vertexID];
    std::shared_ptr<Triangles> trianglesB = meshB->triangles;
    int numTris = trianglesB->triFaceList.size();

    isect->t = std::numeric_limits<T>::infinity();
    isect->hit = false;
    isect->triangleIndex = -1;

    Eigen::Matrix<T, 3, 1> dir = Eigen::Matrix<T, 3, 1>(vCurr - vPrev);
    Ray forwardRay(vPrev, dir); //Ray going towards curr vertex Pos
    float forwardRayLength = forwardRay.direction.norm();
    forwardRay.direction.normalize();
    //-----------------------------------------------------------------------------

    for(int i=0; i < numTris; i++)
    {
        //Test Moving Vertex against Stationary Triangle
        float tTemp = isect->t;
        Eigen::Matrix<T, 3, 1> baryCoords;
        bool intersect = meshB->triangles->intersect(forwardRay, i, &tTemp, meshB->vertices, &baryCoords);

        if((intersect && isect->t > tTemp) && (tTemp <= forwardRayLength))
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp;
            isect->BarycentricWeights = baryCoords;
            isect->normal = meshB->triangles->triNormalList[i];
            isect->point = forwardRay.origin + tTemp * forwardRay.direction;
            isect->penetrationDistance = tTemp * forwardRay.direction;
            isect->vertsOfTriangle = meshB->triangles->triFaceList[i];
            isect->color << 1.0f, 0.0f, 0.0f;
        }

        //------------------------------------------------------//
        //Test Moving Triangle against Stationary Vertex
        if (!isect->hit) 
        {
            Eigen::Matrix<T, 3, 1> triAvgVel;
            trianglesB->computeAvgVelocity(i, triAvgVel, meshB->vertices, meshB->prevVerts);
            Ray rayTVF(vCurr, -triAvgVel); //Ray trying to hit triangle in its previous position with a direction dictated by the -ve of the triangle's velocity
            Ray rayTVB(vCurr, triAvgVel); //Ray trying to hit triangle in its current position with a direction dictated by the triangle's velocity

            //we dont care about the first intersection beyond that it intersects and so we can override 
            //the data in tTemp and baryCoords
            bool intersectF = meshB->triangles->intersect(rayTVF, i, &tTemp, meshB->vertices, &baryCoords);
            bool intersectB = meshB->triangles->intersect(rayTVB, i, &tTemp, meshB->vertices, &baryCoords);
            if(intersectF && intersectB)
            {
                //Point was engulfed by Moving Triangle
                isect->hit = true;
                isect->triangleIndex = i;
                isect->t = tTemp;
                isect->BarycentricWeights = baryCoords;
                isect->normal = meshB->triangles->triNormalList[i];
                isect->point = rayTVB.origin + tTemp * rayTVB.direction;
                isect->penetrationDistance = tTemp * rayTVB.direction;
                isect->vertsOfTriangle = meshB->triangles->triFaceList[i];
                isect->color << 0.0f, 0.0f, 1.0f;
            }
        }
    }

    return isect->hit;
}
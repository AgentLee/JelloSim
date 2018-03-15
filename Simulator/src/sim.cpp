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
    //eulerIntegrationWithCollisionTesting(dt);
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

// returns index of closest triangle on 2nd mesh.. -1 otherwise..
// also sets t to dist of closest tri.. -1 otherwise..
bool Sim::LineTriangleIntersection(const Eigen::Matrix<T, 3, 1>& origPos, const Eigen::Matrix<T, 3, 1>& newPos, std::shared_ptr<Mesh> standingMesh, Intersection *isect, float dt)
{

	Ray r;
	r.origin = origPos;
	r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);
    float length = r.direction.norm();
    r.direction.normalize();

	// assuming only 2 meshes for now..
	int tris = standingMesh->triangles->triFaceList.size();
	isect->t = std::numeric_limits<T>::infinity();
    isect->hit = false;
	isect->triangleIndex = -1;
	for(int i=0; i < tris; i++)
	{
        // CASE 1: POINT IS MOVING THROUGH THE TRIANGLE
		float tTemp = isect->t;
        Eigen::Matrix<T, 3, 1> baryCoords;
		bool intersect = standingMesh->triangles->intersect(r, i, &tTemp, standingMesh->vertices, &baryCoords);

        // CASE 2: TRIANGLE IS MOVING THROUGH THE POINT
        Eigen::Matrix<T, 3, 1> triVel = standingMesh->vertices->vel[standingMesh->triangles->triFaceList[i][0]] +
                                        standingMesh->vertices->vel[standingMesh->triangles->triFaceList[i][1]] +
                                        standingMesh->vertices->vel[standingMesh->triangles->triFaceList[i][2]];
        triVel /= 3.0;
        Ray r1;
        r1.origin = origPos;
        Eigen::Matrix<T, 3, 1> tempPos = origPos + triVel * dt;
        r1.direction = tempPos;
        float length1 = r.direction.norm();
        r1.direction.normalize();
        float tTemp1 = isect->t;
        Eigen::Matrix<T, 3, 1> baryCoords1;
        bool intersect1 = standingMesh->triangles->intersect(r1, i, &tTemp1, standingMesh->vertices, &baryCoords1);

		if(!intersect1 && intersect && isect->t > tTemp && tTemp <= length)
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp;
            isect->BarycentricWeights = baryCoords;
            isect->normal = standingMesh->triangles->triNormalList[i];
            isect->point = r.origin + tTemp * r.direction;
		}
        else if(!intersect && intersect1 && isect->t > tTemp1 && tTemp1 <= length1)
        {
            isect->hit = true;
            isect->triangleIndex = i;
            isect->t = tTemp1;
            isect->BarycentricWeights = baryCoords;
            isect->normal = standingMesh->triangles->triNormalList[i];
            isect->point = r1.origin + tTemp1 * r1.direction;
        }
        else if(intersect && intersect1)
        {
            if(isect->t > tTemp1 && tTemp1 <= length1)
            {
                isect->hit = true;
                isect->triangleIndex = i;
                isect->t = tTemp1;
                isect->BarycentricWeights = baryCoords;
                isect->normal = standingMesh->triangles->triNormalList[i];
                isect->point = r1.origin + tTemp1 * r1.direction;
            }
            else if(isect->t > tTemp && tTemp <= length)
            {
                isect->hit = true;
                isect->triangleIndex = i;
                isect->t = tTemp;
                isect->BarycentricWeights = baryCoords;
                isect->normal = standingMesh->triangles->triNormalList[i];
                isect->point = r.origin + tTemp * r.direction;
            }
        }

	}



//    Ray r1;
//    r1.origin = newPos;
//    r1.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);
//    float length = r.direction.norm();
//    r.direction.normalize();
//
//    // assuming only 2 meshes for now..
//    int tris = standingMesh->triangles->triFaceList.size();
//    for(int i=0; i < tris; i++)
//    {
//        float tTemp = isect->t;
//        Eigen::Matrix<T, 3, 1> baryCoords;
//        bool intersect = standingMesh->triangles->intersect(r, i, &tTemp, standingMesh->vertices, &baryCoords);
//        if(intersect && isect->t > tTemp && tTemp <= length)
//        {
//            isect->hit = true;
//            isect->triangleIndex = i;
//            isect->t = tTemp;
//            isect->BarycentricWeights = baryCoords;
//            isect->normal = MeshList[1]->triangles->triNormalList[i];
//            isect->point = r.origin + tTemp * r.direction;
//        }
//    }

	return isect->hit;
}

//bool Sim::TriangleLineIntersection(const Eigen::Matrix<T, 3, 1>& origPosPoint, int triIdx, const std::shared_ptr<Triangles> triangles, Intersection *isect)
//{
//    // create ray and call triangle intersection for all triangles..
//    // store triangle reference
//    Ray r;
//    r.origin = origPos;
//    r.direction = Eigen::Matrix<T, 3, 1>(newPos - origPos);
//    float length = r.direction.norm();
//    r.direction.normalize();
//
//    // assuming only 2 meshes for now..
//    int tris = MeshList[1]->triangles->triFaceList.size();
//    isect->t = std::numeric_limits<T>::infinity();
//    isect->hit = false;
//    isect->triangleIndex = -1;
//    for(int i=0; i < tris; i++)
//    {
//        float tTemp = isect->t;
//        Eigen::Matrix<T, 3, 1> baryCoords;
//        bool intersect = MeshList[1]->triangles->intersect(r, i, &tTemp, MeshList[1]->vertices, &baryCoords);
//        if(intersect && isect->t > tTemp && tTemp <= length)
//        {
//            isect->hit = true;
//            isect->triangleIndex = i;
//            isect->t = tTemp;
//            isect->BarycentricWeights = baryCoords;
//            isect->normal = MeshList[1]->triangles->triNormalList[i];
//            isect->point = r.origin + tTemp * r.direction;
//        }
//    }
//
//    return isect->hit;
//}

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
        for (uint j =0; j < MeshList.size(); j++)
        {
            if(i==j)
                continue;

            std::shared_ptr<Mesh> fallingMesh = MeshList[i];
            std::shared_ptr<Mesh> standingMesh = MeshList[j];

            bool intersects = Intersect_AABB_with_AABB(fallingMesh->AABB, standingMesh->AABB);

            if (intersects)
            {
                for (int k = 0; k < fallingMesh->vertices->numParticles; ++k)
                {
                    std::shared_ptr<Triangles> triangles = standingMesh->triangles; //i --> current mesh's triangles

                    // for every vertex create a ray from the current position to its projected position in the next frame
                    // See if that ray intersects any triangle that Belongs to the other mesh.
                    //Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[particleIndex] + vertices->vel[particleIndex] * dt;
                    Intersection isect;
                    LineTriangleIntersection(fallingMesh->vertices->pos[k], fallingMesh->prevVerts->pos[k], standingMesh, &isect, dt);

                    if (isect.hit)
                    {
                        //resolve collisions between vertex and a triangle
                        Eigen::Matrix<T, 3, 1> displacement;
                        displacement = fallingMesh->vertices->pos[k] - isect.point;

                        // set position or velocity or both OR use momentum OR use paper's implementation with normal reaction forces and friction
                        Eigen::Matrix<uint, 3, 1> verticesOfTriangle;
                        verticesOfTriangle = triangles->triFaceList[isect.triangleIndex];

                        if (SET_POSITIONS)
                        {
                            //---------------- Setting Positions ---------------------------------
                            //fallingMesh->vertices->pos[k] = fallingMesh->prevVerts->pos[k]; //isect.point.cast<T>() - displacement; // 1/4th to moving vertex
                            newPositions[i][k] = fallingMesh->prevVerts->pos[k] + displacement;
                            //displacement = 0.75f*displacement;				// 3/4th to moving triangle

                            //Eigen::Matrix<T, 3, 1> disp = fallingMesh->vertices->pos[k] - fallingMesh->prevVerts->pos[k];
//                            newPositions[j][verticesOfTriangle[0]] += disp;
//                            newPositions[j][verticesOfTriangle[1]] += disp;
//                            newPositions[j][verticesOfTriangle[2]] += disp;
                        }
                        if (SET_VELOCITIES)
                        {
                            Eigen::Matrix<T, 3, 1> normal;
                            normal = triangles->triNormalList[isect.triangleIndex];

                            Eigen::Matrix<T, 3, 1> d;
                            d = fallingMesh->vertices->vel[k];

                            if(normal.dot(d) < 0.0f)
                            {
                                normal = -normal;
                            }

                            //dir.normalize();
                            Eigen::Matrix<T, 3, 1> reflected;
                            reflected = d - 2.0f*(d.dot(normal))*normal;
                            newVelocities[i][k] = reflected;//Eigen::Matrix<T, 3, 1>::Zero();
//
//                            newVelocities[j][verticesOfTriangle[0]] *= fallingMesh->vertices->vel[k]/3.0;//standingMesh->vertices->vel[verticesOfTriangle[0]] + displacement;
//                            newVelocities[j][verticesOfTriangle[1]] *= fallingMesh->vertices->vel[k]/3.0;
//                            newVelocities[j][verticesOfTriangle[2]] *= fallingMesh->vertices->vel[k]/3.0;
                        }
                    }
                }

//                for (int k = 0; k < fallingMesh->triangles->numTriangles; ++k) {
//                    std::shared_ptr<Triangles> triangles = standingMesh->triangles; //i --> current mesh's triangles
//
//                    // for every vertex create a ray from the current position to its projected position in the next frame
//                    // See if that ray intersects any triangle that Belongs to the other mesh.
//                    //Eigen::Matrix<T, 3, 1> projectedPos = vertices->pos[particleIndex] + vertices->vel[particleIndex] * dt;
//                    Intersection isect;
//                    LineTriangleIntersection(fallingMesh->vertices->pos[k], fallingMesh->prevVerts->pos[k], &isect);
//
//                    if (isect.hit) {
//                        //resolve collisions between vertex and a triangle
//                        Eigen::Matrix<T, 3, 1> displacement;
//                        displacement = fallingMesh->vertices->pos[k] - isect.point;
//
//                        // set position or velocity or both OR use momentum OR use paper's implementation with normal reaction forces and friction
//                        Eigen::Matrix<uint, 3, 1> verticesOfTriangle;
//                        verticesOfTriangle = triangles->triFaceList[isect.triangleIndex];
//
//                        if (SET_POSITIONS)
//                        {
//                            //---------------- Setting Positions ---------------------------------
//                            //fallingMesh->vertices->pos[k] = fallingMesh->prevVerts->pos[k]; //isect.point.cast<T>() - displacement; // 1/4th to moving vertex
//                            newPositions[i][k] = fallingMesh->prevVerts->pos[k] - 1.3 * displacement;
//                            //displacement = 0.75f*displacement;				// 3/4th to moving triangle
//
//                            // move vertices of triangle according to barycentric weights
////                            newPositions[j][verticesOfTriangle[0]] = standingMesh->vertices->pos[verticesOfTriangle[0]] + displacement;
////                            newPositions[j][verticesOfTriangle[1]] = standingMesh->vertices->pos[verticesOfTriangle[1]] + displacement;
////                            newPositions[j][verticesOfTriangle[2]] = standingMesh->vertices->pos[verticesOfTriangle[2]] + displacement;
//                        }
//                        if (SET_VELOCITIES)
//                        {
//                            Eigen::Matrix<T, 3, 1> normal;
//                            normal = triangles->triNormalList[isect.triangleIndex];
//
//                            Eigen::Matrix<T, 3, 1> d;
//                            d = fallingMesh->vertices->vel[k];
//
//                            if(normal.dot(d) < 0.0f)
//                            {
//                                normal = -normal;
//                            }
//
//                            //dir.normalize();
//                            Eigen::Matrix<T, 3, 1> reflected;
//                            reflected = d - 2.0f*(d.dot(normal))*normal;
//                            newVelocities[i][k] = reflected;//Eigen::Matrix<T, 3, 1>::Zero();
//
////                            newVelocities[j][verticesOfTriangle[0]] = standingMesh->vertices->vel[verticesOfTriangle[0]] + fallingMesh->vertices->vel[k];//standingMesh->vertices->vel[verticesOfTriangle[0]] + displacement;
////                            newVelocities[j][verticesOfTriangle[1]] = standingMesh->vertices->vel[verticesOfTriangle[1]] + fallingMesh->vertices->vel[k];
////                            newVelocities[j][verticesOfTriangle[2]] = standingMesh->vertices->vel[verticesOfTriangle[2]] + fallingMesh->vertices->vel[k];
//                        }
//                    }
//                }
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

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>

#include "particles.h"
#include "triangles.h"
#include "tetrahedrons.h"
#include "bounds.h"

class Mesh
{
public:
    std::shared_ptr<Particles> prevVerts;
	std::shared_ptr<Particles> vertices;
	std::shared_ptr<Triangles> triangles;
	std::shared_ptr<Tetrahedrons> tetras;
	Bounds AABB;

	// grid that stores the triangles of the mesh; 
	// grid extends beyond the biggest AABB that could be formed to accommodate stretching and
	// prevent resizing of the grid
	std::vector<uint> grid;

	Mesh( const std::string& nodeFile, const std::string& faceFile, 
		const std::string& eleFile, const std::string& objFile,
		float _gridCellSize, float _density, float _poissonsRatio, float _youngsModulus );
	Mesh( std::shared_ptr<Particles> _vertices, std::shared_ptr<Triangles> _triangles, 
		std::shared_ptr<Tetrahedrons> _tetras, float _gridCellSize );

	// Method of computing a 1D index from a 3D grid index.
	int gridIndex3Dto1D(int x, int y, int z, int gridResolution);

	void calcBounds();
	void updateBounds(const Bounds &newAABB);
	void initializeGrid();

	void initMeshForSim();
	void computeElasticForcesOnMesh(int frame);
	void clearForces();

	void resetMass();
	void translateMesh(Vector3f& translation);

    void setPrevVerts();

	void tetRead( const std::string& nodeFile, const std::string& faceFile, const std::string& eleFile, const std::string& objFile );

private:
	float gridResolution;
};
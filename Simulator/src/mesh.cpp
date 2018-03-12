#include "mesh.h"

Mesh::Mesh( const std::string& nodeFile, const std::string& faceFile, 
			const std::string& eleFile, const std::string& objFile,
			float _gridCellSize, float _density, float _poissonsRatio, float _youngsModulus ) 
		: gridResolution(_gridCellSize)
{
	vertices  = std::make_shared<Particles>(0, 0.0f);
	triangles = std::make_shared<Triangles>();
	tetras    = std::make_shared<Tetrahedrons>( _density, _poissonsRatio, _youngsModulus );

	tetRead( nodeFile, faceFile, eleFile, objFile );
	resetMass();

    prevVerts  = std::make_shared<Particles>(vertices->numParticles, 0.0f);

	AABB = Bounds();
	calcBounds();
	initializeGrid();
}

Mesh::Mesh( std::shared_ptr<Particles> _vertices, std::shared_ptr<Triangles> _triangles, 
			std::shared_ptr<Tetrahedrons> _tetras, float _gridCellSize )
		: vertices(_vertices), triangles(_triangles), tetras(_tetras), gridResolution(_gridCellSize)
{
	resetMass();

	AABB = Bounds();
	calcBounds();
	initializeGrid();
}

// Method of computing a 1D index from a 3D grid index.
int Mesh::gridIndex3Dto1D(int x, int y, int z, int gridResolution)
{
	return x + y * gridResolution + z * gridResolution * gridResolution;
}

void Mesh::resetMass()
{
	for(int i=0; i<vertices->numParticles; i++)
	{
		vertices->mass[i] = 0.0f;
	}
}

void Mesh::initMeshForSim()
{
	// Precompute rest deformation (Dm), volume, inverse Dm, and volume*inverseDmTranspose for each tetrahedron
	tetras->restDeformation.resize(tetras->numTetra);
	tetras->restInverseDeformation.resize(tetras->numTetra);
	tetras->undeformedVolume.resize(tetras->numTetra);
	tetras->undefVol_into_restInvDefTranspose.resize(tetras->numTetra);

	for(int j=0; j<tetras->numTetra; j++)
	{
		tetras->computeRestDeformation( j, vertices );
		tetras->computeInvRestDeformation( j );
		tetras->computeUndeformedVolume( j );
		tetras->computeUndefVol_into_restInvDefTranspose( j );

		tetras->addMass( j, vertices );
	}
}

void Mesh::computeElasticForcesOnMesh(int frame)
{
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

void Mesh::clearForces()
{
	//Set forces for all vertices/particles to zero
	for( int j=0; j < vertices->numParticles; j++ )
	{
		vertices->force[j] = Eigen::Matrix<T, 3, 1>::Zero();
	}
}

//separate function because we have to tell the bounding box to move too
void Mesh::translateMesh(Vector3f& translation)
{
	for(int i=0; i<vertices->numParticles; i++)
	{
		vertices->pos[i] += translation.cast<T>();
	}

	AABB.min += translation;
	AABB.max += translation;
}

void Mesh::calcBounds()
{
	//loop over vertices to identify the Bounding Box
    double inf = std::numeric_limits<T>::infinity();
    AABB.min << inf, inf, inf;
    AABB.max << -inf, -inf, -inf;

	for(int i=0; i<vertices->numParticles; i++)
	{
		Point3f pos = vertices->pos[i].cast<float>();
		AABB = Union( AABB, pos );
	}
}

void Mesh::updateBounds(const Bounds &newAABB)
{
	//AABB.min = std
}

void Mesh::initializeGrid()
{

}

// Reads all the required files and stores stuff in respective arrays
void Mesh::tetRead( const std::string& nodeFile, const std::string& faceFile, const std::string& eleFile, const std::string& objFile )
{
	// read node file and store list of points (id, pos, attributes)
	// read .ele file and store tetrahedra (id, nodes)
	// can read more files for faces, edges, etc.

	vertices->tetgen_readNode( nodeFile ); //read vertices
	triangles->tetgen_readFace( faceFile );
	tetras->tetgen_readEleFile( eleFile ); //read in tetrahedrons with particle indices
	triangles->create_objFile( objFile, vertices );
}

void Mesh::setPrevVerts()
{
    for(int i=0; i<prevVerts->numParticles; i++)
    {
        prevVerts->pos[i] = vertices->pos[i];
        prevVerts->vel[i] = vertices->vel[i];
    }
}


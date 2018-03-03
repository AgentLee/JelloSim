#include <Partio.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <string>

#include <memory>
#include <iostream>
#include <time.h>

#include "src/particles.h"
#include "src/triangles.h"
#include "src/tetrahedrons.h"
#include "src/mesh.h"
#include "src/utilities.h"
#include "src/sim.h"

using T = double;

constexpr int beginFrame = 2; //1 is initial state that is written separately
constexpr int endFrame = 30;
constexpr int numSimulationStepsPerFrame = 400;
constexpr float dt = 1e-4;

const std::string cube_nodeFile = "../Assets/Meshes/cube_poly_0.001/cube.1.node";
const std::string cube_faceFile = "../Assets/Meshes/cube_poly_0.001/cube.1.face";
const std::string cube_eleFile  = "../Assets/Meshes/cube_poly_0.001/cube.1.ele";
const std::string cube_objFile  = "../Assets/OBJs/cube1.obj";
const std::string baseFileNamePoints = "../Assets/BGEOs/jelloTestFrame";

void createScene( std::vector<std::shared_ptr<Mesh>>& MeshList )
{
	const float gridCellSize = 1.0f;

	std::shared_ptr<Mesh> cube = std::make_shared<Mesh>( gridCellSize );
	Utils::tetRead( cube->vertices, cube->triangles, cube->tetras, 
					cube_nodeFile, cube_faceFile, cube_eleFile, cube_objFile );

    for(int i=0; i<cube->vertices->numParticles; i++)
    {
    	//move mesh up by 1m and reset mass to zero
        cube->vertices->pos[i](1) += 1.0f;
        cube->vertices->mass[i] = 0.0f;
    }

    MeshList.push_back(cube);
}

int main(int argc, char* argv[])
{
	std::vector<std::shared_ptr<Mesh>> MeshList;
	createScene( MeshList );

	//Create Bgeo file for first frame
	std::cout << "Simualting Frame: 1" << "..." << std::endl;
	Utils::writeBGEOforFrame( baseFileNamePoints, 1, MeshList );

	// Start sim
	//Initialize Jello Sim
	clock_t t;
	t = clock();

	std::shared_ptr<Sim> sim = std::make_shared<Sim>( MeshList[0]->tetras, MeshList[0]->vertices );
	sim->init();

	bool collided = false;

	//Main Loop of Jello Sim
	for(int i=beginFrame; i<=endFrame; i++)
	{
		std::cout << "Simualting Frame: " << i << "..." << std::endl;
		for( int j = 0; j <= numSimulationStepsPerFrame; j++ )
		{
			sim->update(dt, i, collided);
		}

		Utils::writeBGEOforFrame( baseFileNamePoints, i, MeshList );
	}
	t = clock() - t;	// End sim

	std::cout << "Simulated in " << ((float)t)/CLOCKS_PER_SEC << " seconds." << std::endl;
	
	return 0;
}

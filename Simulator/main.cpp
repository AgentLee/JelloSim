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
constexpr int endFrame = 70;
constexpr int numSimulationStepsPerFrame = 300;
constexpr float dt = 1e-4;
constexpr float density = 1000.f;
constexpr float youngsModulus = 500000.0f;
constexpr float poissonsRatio = 0.3f;

const std::string nodeFileNames[] = { "../Assets/Meshes/cube_poly_0.001/cube.1.node",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.node" };
const std::string faceFileNames[] = { "../Assets/Meshes/cube_poly_0.001/cube.1.face",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.face" };
const std::string eleFileNames[]  = { "../Assets/Meshes/cube_poly_0.001/cube.1.ele",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.ele"};
const std::string objFileNames[]  = { "../Assets/OBJs/FirstCube.obj",
									  "../Assets/OBJs/SecondCube.obj" };
const std::string bgeoFileNames[] = { "../Assets/BGEOs/jelloCube1Frame",
									  "../Assets/BGEOs/jelloCube2Frame" };

void createScene( std::vector<std::shared_ptr<Mesh>>& MeshList )
{
	const float gridCellSize = 1.0f;
	Vector3f translation = Vector3f(0.0f, 1.0f, 0.0f);

	{
		std::shared_ptr<Mesh> cube1 = std::make_shared<Mesh>( nodeFileNames[0], faceFileNames[0], 
															eleFileNames[0], objFileNames[0], gridCellSize, 
															density, poissonsRatio, youngsModulus );
		cube1->translateMesh(translation);
		MeshList.push_back(cube1);
	}

	{
		translation = Vector3f(3.0f, 3.0f, 0.0f);
		std::shared_ptr<Mesh> cube2 = std::make_shared<Mesh>( nodeFileNames[1], faceFileNames[1], 
															eleFileNames[1], objFileNames[1], gridCellSize, 
															density, poissonsRatio, youngsModulus );
		cube2->translateMesh(translation);
		MeshList.push_back(cube2);
	}
}

void writeBGEOsforMeshes( std::vector<std::shared_ptr<Mesh>>& MeshList, int frameNumber )
{
	std::cout << "Simualting Frame: " << frameNumber << "..." << std::endl;
	for(uint i=0; i<MeshList.size(); i++)
	{
		Utils::writeBGEOforFrame( bgeoFileNames[i], frameNumber, MeshList[i]->vertices );
	}
}

int main(int argc, char* argv[])
{
	std::vector<std::shared_ptr<Mesh>> MeshList;
	createScene( MeshList );

	//Create Bgeo file for first frame
	writeBGEOsforMeshes( MeshList, 1 );

	// Start sim
	// Initialize Jello Sim
	clock_t t;
	t = clock();
	
	std::shared_ptr<Sim> sim = std::make_shared<Sim>( MeshList );
	sim->init();

	bool collided = false;

	std::cout << "Num Meshes " << MeshList.size() << std::endl;

	//Main Loop of Jello Sim
	for(int i=beginFrame; i<=endFrame; i++)
	{
		for( int j = 0; j <= numSimulationStepsPerFrame; j++ )
		{
			sim->update(dt, i, collided);
		}

		writeBGEOsforMeshes( MeshList, i );
	}
	t = clock() - t;	// End sim

	std::cout << "Simulated in " << ((float)t)/CLOCKS_PER_SEC << " seconds." << std::endl;
	
	return 0;
}

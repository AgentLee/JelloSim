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

//#define CONVERT_OBJ_TO_POLY

using T = double;

constexpr int beginFrame = 2; //1 is initial state that is written separately
constexpr int endFrame = 120;
constexpr int numSimulationStepsPerFrame = 300;
constexpr float dt = 1e-4;
constexpr float density = 1000.f;
constexpr float youngsModulus = 500000.0f;
constexpr float poissonsRatio = 0.3f;

const std::string nodeFileNames[] = { "../Assets/Meshes/cube_poly_0.0625/cube.1.node",
									  "../Assets/Meshes/cube_poly_0.0625/cube.1.node" };
const std::string faceFileNames[] = { "../Assets/Meshes/cube_poly_0.0625/cube.1.face",
									  "../Assets/Meshes/cube_poly_0.0625/cube.1.face" };
const std::string eleFileNames[]  = { "../Assets/Meshes/cube_poly_0.0625/cube.1.ele",
									  "../Assets/Meshes/cube_poly_0.0625/cube.1.ele" };
const std::string objFileNames[]  = { "../Assets/OBJs/FirstCube.obj",
									  "../Assets/OBJs/SecondCube.obj" };
const std::string bgeoFileNames[] = { "../Assets/BGEOs/jelloCube1Frame",
									  "../Assets/BGEOs/jelloCube2Frame" };

const std::string objToPolyNames[] = { "../Assets/objs_polys/teapotURN.obj",
                                       "../Assets/objs_polys/teapotURN.poly" };

void createScene( std::vector<std::shared_ptr<Mesh>>& MeshList, Bounds* FixedRegion )
{
	const float gridCellSize = 1.0f;
	Vector3f translation = Vector3f(0.1f, 1.6f, 0.1f);

	{
		std::shared_ptr<Mesh> cube1 = std::make_shared<Mesh>( nodeFileNames[0], faceFileNames[0], 
															eleFileNames[0], objFileNames[0], gridCellSize, 
															density, poissonsRatio, youngsModulus );
		cube1->translateMesh(translation);
		MeshList.push_back(cube1);
	}

	{
		translation = Vector3f(0.0f, 0.0f, 0.0f);
		std::shared_ptr<Mesh> cube2 = std::make_shared<Mesh>( nodeFileNames[1], faceFileNames[1], 
															eleFileNames[1], objFileNames[1], gridCellSize, 
															density, poissonsRatio, youngsModulus );
		cube2->translateMesh(translation);
		MeshList.push_back(cube2);
	}

    // // Fixed point region
    // {
    //     FixedRegion->min << -0.5, 8.0, -1.0; // good values for bunneh
    //     FixedRegion->max << 0.0, 9.0, 1.0;
    // }
}

void writeBGEOsforMeshes( std::vector<std::shared_ptr<Mesh>>& MeshList, int frameNumber )
{
	for(uint i=0; i<MeshList.size(); i++)
	{
		Utils::writeBGEOforFrame( bgeoFileNames[i], frameNumber, MeshList[i]->vertices );
	}
	std::cout << "Frame " << frameNumber << " Simualted "  << "..." << std::endl;
}

int main(int argc, char* argv[])
{
#ifdef CONVERT_OBJ_TO_POLY
    Utils::objToPoly(objToPolyNames[0], objToPolyNames[1]);
#else
	std::vector<std::shared_ptr<Mesh>> MeshList;
    Bounds FixedRegion;
	createScene( MeshList, &FixedRegion );

	//Create Bgeo file for first frame
	writeBGEOsforMeshes( MeshList, 1 );

	// Start sim
	// Initialize Jello Sim
	clock_t t;
	t = clock();
	
	std::shared_ptr<Sim> sim = std::make_shared<Sim>( MeshList, FixedRegion );
	sim->init();

	std::cout << "Num Meshes " << MeshList.size() << std::endl;

	//Main Loop of Jello Sim
	for(int i=beginFrame; i<=endFrame; i++)
	{
		for( int j = 0; j <= numSimulationStepsPerFrame; j++ )
		{
			sim->update(dt, i);
		}

		writeBGEOsforMeshes( MeshList, i );
	}
	t = clock() - t;	// End sim

	std::cout << "Simulated in " << ((float)t)/CLOCKS_PER_SEC << " seconds." << std::endl;
#endif
	return 0;
}

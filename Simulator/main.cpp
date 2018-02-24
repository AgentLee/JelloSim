#include <Partio.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <string>

#include <memory>
#include <iostream>

#include "src/particles.h"
#include "src/triangles.h"
#include "src/tetrahedrons.h"
#include "src/utilities.h"
#include "src/sim.h"

using T = double;

constexpr int beginFrame = 2; //1 is initial state that is written separately
constexpr int endFrame = 120;
constexpr int numSimulationStepsPerFrame = 10;

const std::string nodeFile = "../Assets/Meshes/cube_poly_0.001/cube.1.node";
const std::string faceFile = "../Assets/Meshes/cube_poly_0.001/cube.1.face";
const std::string eleFile  = "../Assets/Meshes/cube_poly_0.001/cube.1.ele";
const std::string objFile  = "../Assets/OBJs/cube1.obj";
const std::string baseFileNamePoints = "../Assets/BGEOs/jelloTestFrame";

int main(int argc, char* argv[])
{
	std::cout << "BEGIN SIM" << std::endl;

	std::shared_ptr<Particles> vertices = std::make_shared<Particles>(0, 0.1f);
	std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();
	std::shared_ptr<Tetrahedrons> tetras = std::make_shared<Tetrahedrons>();

	Utils::tetRead( vertices, triangles, tetras, nodeFile, faceFile, eleFile, objFile );

    for(int i=0; i<vertices->numParticles; i++)
    {
        vertices->pos[i](1) += .5; //computes and adds elastic forces to each particle
    }

	//Create Bgeo file for first frame
	Utils::writeBGEOforFrame( baseFileNamePoints, 1, vertices );

	//Initialize Jello Sim
	std::shared_ptr<Sim> sim = std::make_shared<Sim>( tetras, vertices );
	sim->init();

	bool collided = false;

	//Main Loop of Jello Sim
	for(int i=beginFrame; i<=endFrame; i++)
	{
		for(int j=0; j<=numSimulationStepsPerFrame; j++)
		{
			sim->update(0.00001, i, collided);
		}

		Utils::writeBGEOforFrame( baseFileNamePoints, i, vertices );
	}

	std::cout << "END SIM" << std::endl;
	return 0;
}

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
constexpr int dim = 3;
constexpr int beginFrame = 2; //1 is initial state that is written separately
constexpr int endFrame = 120;

int main(int argc, char* argv[])
{
	std::shared_ptr<Particles> vertices = std::make_shared<Particles>(0, 1.0f);
	std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();
	std::shared_ptr<Tetrahedrons> tetras = std::make_shared<Tetrahedrons>();

	//read vertices
	std::string fileName = "../Meshes/cube_poly_1/cube.1.node";
	vertices->tetgen_readNode( fileName );

	//Create Bgeo file for first frame
	std::string baseFileNamePoints = "TestPoints/testPointsFrame";
	std::string pointsFile = baseFileNamePoints + "1.bgeo";
	Utils::writePartio<T,dim>(pointsFile, vertices, vertices->numParticles);

	//read faces and create obj file
	std::string faceFile = "../Meshes/cube_poly_1/cube.1.face";
	std::string objFileName = "../Meshes/cube_poly_1/cube.1.node";
	Utils::generateTriObjFile( triangles, vertices, fileName, objFileName );

	//read in tetrahedrons with particle indices
	std::string eleFile = "../Meshes/cube_poly_1/cube.1.ele";
	tetras->tetgen_readEleFile(eleFile);
	
	//Initialize Jello Sim
	std::shared_ptr<Sim> sim = std::make_shared<Sim>( tetras, vertices );
	sim->init();

	for(int i=beginFrame; i<=endFrame; i++)
	{
		sim->update();
		//create bgeo file for current frame
	}

	return 0;
}

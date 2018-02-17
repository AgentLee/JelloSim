#include <Partio.h>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <memory>
#include <iostream>

#include "src/particles.h"
#include "src/triangles.h"

#include "src/utilities.h"
#include "src/forces.h"

using T = double;
constexpr int dim = 3;
constexpr int numVertices = 100;
constexpr int numFrames = 120;

int main(int argc, char* argv[])
{
	std::shared_ptr<Particles> vertices = std::make_shared<Particles>(0, 1.0f);
	std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();

	std::string fileName = "../Meshes/cube_poly_1/cube.1.node";
	vertices->tetgen_readNode( fileName );

	std::string baseFileNamePoints = "TestPoints/testPointsFrame";
	std::string pointsFile = baseFileNamePoints + "1.bgeo";
	Utils::writePartio<T,dim>(pointsFile, vertices, vertices->numParticles);

	std::string faceFile = "../Meshes/cube_poly_1/cube.1.face";
	std::string objFileName = "../Meshes/cube_poly_1/cube.1.node";
	Utils::generateTriObjFile( triangles, vertices, fileName, objFileName );

	return 0;
}

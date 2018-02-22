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
const std::string baseFileNamePoints = "TestPoints/testPointsFrame";

const float timeStepsPerFrame = 10.0;
const float frameRate = 24.0;

int main(int argc, char* argv[])
{
	std::shared_ptr<Particles> vertices = std::make_shared<Particles>(0, 1.0f);
	std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();
	std::shared_ptr<Tetrahedrons> tetras = std::make_shared<Tetrahedrons>();

	Utils::tetRead( vertices, triangles, tetras );

    for(int i=0; i<vertices->numParticles; i++)
    {
        vertices->pos[i](1) += 5.0; //computes and adds elastic forces to each particle
    }

	//Create Bgeo file for first frame
	//Utils::writeBGEOforFrame( baseFileNamePoints, 1, vertices );
	std::string pointsFile = baseFileNamePoints;
	pointsFile += std::to_string(1);
	pointsFile += ".bgeo";
	Utils::writePartio<T, 3>(pointsFile, vertices, vertices->numParticles);

	//Initialize Jello Sim
	std::shared_ptr<Sim> sim = std::make_shared<Sim>( tetras, vertices );
	sim->init();

	//Main Loop of Jello Sim
	for(int i=beginFrame; i<=endFrame; i++)
	{
		for(int j=0; j<=timeStepsPerFrame; j++)
		{
			sim->update(1 / frameRate / timeStepsPerFrame);
		}
		//Utils::writeBGEOforFrame( baseFileNamePoints, i, vertices );
		std::string pointsFile = baseFileNamePoints;
		pointsFile += std::to_string(i);
		pointsFile += ".bgeo";
		Utils::writePartio<T, 3>(pointsFile, vertices, vertices->numParticles);
	}

	return 0;
}

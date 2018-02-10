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

void generatePointsFromParticles( T deltaTime, std::shared_ptr<Particles> particles, std::string& fileName )
{
	//generate points
	std::string index;
	std::string file;
	for(int i=1; i<numFrames; i++)
	{
		index = std::to_string(i);
		file = fileName + index + ".bgeo";

		//update vertices
		particles->updateParticlePositions(deltaTime);
		Utils::writePartio<T,dim>(file, particles, numVertices);
	}
}

void testGenerationOfPoints_LineSegments_and_Triangles()
{
	std::shared_ptr<Particles> particles = std::make_shared<Particles>(numVertices, 1.0f);
	particles->initializeParticles_random();

	std::string baseFileNamePoints = "TestPoints/testPointsFrame";
	std::string testPointsFile = baseFileNamePoints + "0.bgeo";
	Utils::writePartio<T,dim>(testPointsFile, particles, numVertices);
	T deltaTime = 0.1f;

	generatePointsFromParticles( deltaTime, particles, baseFileNamePoints );
	std::string triangleObjFile = "trianglesData.obj";
	Utils::generateTrianglesFromParticles( deltaTime, particles, triangleObjFile );
}

std::vector<T> readSpringForcesInputFile( std::string& file, int& numTestCases )
{
	std::vector<T> inputFloats;
	string line; //cause the input text file is formatted as a giant line made of floats, it is all in a single line
	ifstream springForcesInputFile (file);

	if (springForcesInputFile.is_open())
	{
		while ( getline (springForcesInputFile, line) )
		{
			std::vector<T> lineOfFloatsAsVector = Utils::split<T>(line);
			inputFloats.insert(std::end(inputFloats), std::begin(lineOfFloatsAsVector), std::end(lineOfFloatsAsVector));

			numTestCases++;
		}
		springForcesInputFile.close();
	}

	return inputFloats;
}

void writeSpringForcesOutputFile( ofstream& outputFile, std::vector<Eigen::Matrix<T,3,1>> outForces )
{
	for(uint i=0; i<outForces.size(); i++)
	{
		outputFile << outForces[i][0] << " " << outForces[i][1] << " " << outForces[i][2] << " ";
	}

	outputFile << "\n";
}

void checkTestCase( ofstream& outFile, T youngsModulus, T dampingCoeff, T springRestLength, std::shared_ptr<Particles> particles )
{
	//set positions and velocities for the Points -- Test Case 1
	std::vector<Eigen::Matrix<T,3,1>> outForces;

	std::vector<Eigen::Matrix<T,3,1>> springForces;
	std::vector<Eigen::Matrix<T,3,1>> dampingForces;
	//Solve Spring Forces and Damping Forces
	springForces.push_back( Forces::getSpringForce( youngsModulus, springRestLength, particles->pos[0], particles->pos[1] ) );
	springForces.push_back( Forces::getSpringForce( youngsModulus, springRestLength, particles->pos[1], particles->pos[2] ) );
	springForces.push_back( Forces::getSpringForce( youngsModulus, springRestLength, particles->pos[2], particles->pos[3] ) );
	springForces.push_back( Forces::getSpringForce( youngsModulus, springRestLength, particles->pos[3], particles->pos[0] ) );
	springForces.push_back( Forces::getSpringForce( youngsModulus, springRestLength, particles->pos[1], particles->pos[3] ) );

	outForces.push_back(  springForces[0] - springForces[3] );
	outForces.push_back( -springForces[0] + springForces[1] + springForces[4] );
	outForces.push_back( -springForces[1] + springForces[2] );
	outForces.push_back( -springForces[2] + springForces[3] - springForces[4] );

	dampingForces.push_back( Forces::getSpringDampingForce( dampingCoeff, particles->pos[0], particles->pos[1], particles->vel[0], particles->vel[1] ) );
	dampingForces.push_back( Forces::getSpringDampingForce( dampingCoeff, particles->pos[1], particles->pos[2], particles->vel[1], particles->vel[2] ) );
	dampingForces.push_back( Forces::getSpringDampingForce( dampingCoeff, particles->pos[2], particles->pos[3], particles->vel[2], particles->vel[3] ) );
	dampingForces.push_back( Forces::getSpringDampingForce( dampingCoeff, particles->pos[3], particles->pos[0], particles->vel[3], particles->vel[0] ) );
	dampingForces.push_back( Forces::getSpringDampingForce( dampingCoeff, particles->pos[1], particles->pos[3], particles->vel[1], particles->vel[3] ) );

	outForces.push_back( dampingForces[0] - dampingForces[3] );
	outForces.push_back( -dampingForces[0] + dampingForces[1] + dampingForces[4] );
	outForces.push_back( -dampingForces[1] + dampingForces[2] );
	outForces.push_back( -dampingForces[2] + dampingForces[3] - dampingForces[4] );

	//Store Forces in output file
	writeSpringForcesOutputFile( outFile, outForces );
}

int main(int argc, char* argv[])
{
	int numTestCases = 0;
	std::string file = "InputOutput/actual_input.txt";
	std::vector<T> inputFloats = readSpringForcesInputFile(file, numTestCases);

	file = "InputOutput/actual_output.txt";
	ofstream springForcesOutputFile (file);

	int numPoints = 4;
	if (springForcesOutputFile.is_open())
	{
		int count = 0;
		for( int j=0; j<numTestCases; j++ )
		{
			std::shared_ptr<Particles> particles = std::make_shared<Particles>(numPoints, 1.0f);

			for(int i=0; i<numPoints; i++)
			{
				particles->mass[i] = 1;
				particles->pos[i] = ( Eigen::Matrix<T,3,1>() << inputFloats[i*3+3 + count], 
																inputFloats[i*3+4 + count], 
																inputFloats[i*3+5 + count] ).finished();
				particles->vel[i] = ( Eigen::Matrix<T,3,1>() << inputFloats[i*3+3+3*numPoints + count], 
																inputFloats[i*3+4+3*numPoints + count], 
																inputFloats[i*3+5+3*numPoints + count] ).finished();
			}
			
			T youngsModulus = inputFloats[0 + count];
			T dampingCoeff = inputFloats[1 + count];
			T springRestLength = inputFloats[2 + count];

			checkTestCase( springForcesOutputFile, youngsModulus, dampingCoeff, springRestLength, particles );

			count += 27;
		}

		springForcesOutputFile.close();
	}

	return 0;
}

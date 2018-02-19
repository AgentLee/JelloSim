#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdlib.h>
#include <time.h>

#include "particles.h"
#include "triangles.h"

using namespace std;
using T = double;

namespace Utils
{
	template <class T, int dim>
	void inline writePartio(std::string& particleFile, std::shared_ptr<Particles> vertices, int numVertices)
	{
		Partio::ParticlesDataMutable* parts = Partio::create();
		Partio::ParticleAttribute posH, vH, mH;
		mH = parts->addAttribute("m", Partio::VECTOR, 1);
		posH = parts->addAttribute("position", Partio::VECTOR, 3);
		vH = parts->addAttribute("v", Partio::VECTOR, 3);

		for (int i=0; i<numVertices; i++)
		{
			int idx = parts->addParticle();
			float* m = parts->dataWrite<float>(mH, idx);
			float* p = parts->dataWrite<float>(posH, idx);
			float* v = parts->dataWrite<float>(vH, idx);

			m[0] = vertices->mass[0];
			for (uint k = 0; k < 3; k++)
			{
				p[k] = vertices->pos[i][k];
				v[k] = vertices->vel[i][k];
			}
		}

		Partio::write(particleFile.c_str(), *parts);
		parts->release();
	};

	void inline create_objFile(std::string file_name, std::shared_ptr<Particles> vertices, std::shared_ptr<Triangles> triangles)
	{
		ofstream objfile;
		objfile.open (file_name);
		objfile << "# " + file_name + "\n"; //comment with file name

		objfile << "# vertices\n";
		for(int i=0; i<vertices->numParticles; i++)
		{
			objfile << "v  ";
			for (uint k = 0; k < 3; k++)
			{
				objfile << std::to_string(vertices->pos[i][k]) + "  ";
			}
			objfile << "\n";
		}

		objfile << "# faces\n";
		for(int i=0; i<triangles->numTriangles; i++)
		{
			objfile << "f  " + std::to_string(triangles->triFaceList[i][0] + 1) + "//  "
							 + std::to_string(triangles->triFaceList[i][1] + 1) + "//  "
							 + std::to_string(triangles->triFaceList[i][2] + 1) + "//  " + "\n";
		}

		objfile << "# end of obj file\n";
		objfile.close();
	}

	void inline generateTrianglesFromParticles( T deltaTime, std::shared_ptr<Particles> particles, std::string fileName )
	{
		//generate obj file
		std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();

		for(int i=0; i<particles->numParticles-2; i+=2)
		{
			triangles->addTriangle(i, i+1, i+2);
		}
		Utils::create_objFile(fileName, particles, triangles);
	}

	void inline generateTriObjFile( std::shared_ptr<Triangles> triangles,  std::shared_ptr<Particles> particles, const std::string& fileName, const std::string& objFileName )
	{
		//generate obj file
		triangles->tetgen_readFace(fileName);
		Utils::create_objFile(objFileName, particles, triangles);
	}

	// Reads all the required files and stores stuff in respective arrays
	void inline tetRead( std::shared_ptr<Particles> vertices, std::shared_ptr<Triangles> triangles, std::shared_ptr<Tetrahedrons> tetras )
	{
		// read node file and store list of points (id, pos, attributes)
		// read .ele file and store tetrahedra (id, nodes)
		// can read more files for faces, edges, etc.

		//read vertices
		std::string fileName = "../Meshes/cube_poly_1/cube.1.node";
		vertices->tetgen_readNode( fileName );

		//read faces and create obj file
		std::string faceFile = "../Meshes/cube_poly_1/cube.1.face";
		std::string objFileName = "../Meshes/cube_poly_1/cube.1.node";
		Utils::generateTriObjFile( triangles, vertices, fileName, objFileName );

		//read in tetrahedrons with particle indices
		std::string eleFile = "../Meshes/cube_poly_1/cube.1.ele";
		tetras->tetgen_readEleFile(eleFile);
	}

	void inline writeBGEOforFrame( std::string baseFileNamePoints, uint frameNumber, std::shared_ptr<Particles> vertices )
	{
		//create bgeo file for current frame
		std::string pointsFile = baseFileNamePoints + std::to_string(frameNumber) + ".bgeo";
		Utils::writePartio<T, 3>(pointsFile, vertices, vertices->numParticles);
	}

	template<typename T>
	std::vector<T> split(const std::string& line) 
	{
		std::istringstream is(line);
		return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
	}
};
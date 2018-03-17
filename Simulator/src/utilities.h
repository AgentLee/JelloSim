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
		Partio::ParticleAttribute posH, vH, mH, col;
		mH = parts->addAttribute("m", Partio::VECTOR, 1);
		posH = parts->addAttribute("position", Partio::VECTOR, 3);
		vH = parts->addAttribute("v", Partio::VECTOR, 3);
		col = parts->addAttribute("c", Partio::VECTOR, 3);

		for (int i=0; i<numVertices; i++)
		{
			int idx = parts->addParticle();
			float* m = parts->dataWrite<float>(mH, idx);
			float* p = parts->dataWrite<float>(posH, idx);
			float* v = parts->dataWrite<float>(vH, idx);
			float* c = parts->dataWrite<float>(col, idx);

			m[0] = vertices->mass[0];
			for (uint k = 0; k < 3; k++)
			{
				p[k] = vertices->pos[i][k];
				v[k] = vertices->vel[i][k];
				c[k] = vertices->color[i][k];
			}
		}

		Partio::write(particleFile.c_str(), *parts);
		parts->release();
	};

	void inline create_objFile(std::string file_name, std::shared_ptr<Particles>& vertices, std::shared_ptr<Triangles>& triangles)
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

	void inline generateTrianglesFromParticles( T& deltaTime, std::shared_ptr<Particles>& particles, std::string& fileName )
	{
		//generate obj file
		std::shared_ptr<Triangles> triangles = std::make_shared<Triangles>();

		for(int i=0; i<particles->numParticles-2; i+=2)
		{
			triangles->addTriangle(i, i+1, i+2);
		}
		Utils::create_objFile(fileName, particles, triangles);
	}

	void inline generateTriObjFile( std::shared_ptr<Triangles>& triangles,  std::shared_ptr<Particles>& particles, const std::string& fileName, const std::string& objFileName )
	{
		//generate obj file
		triangles->tetgen_readFace(fileName);
		Utils::create_objFile(objFileName, particles, triangles);
	}

	void inline writeBGEOforFrame( const std::string bgeoFileName, uint frameNumber, std::shared_ptr<Particles>& vertices )
	{
		//create bgeo file for current frame
		std::string pointsFile = bgeoFileName;
		pointsFile += std::to_string(frameNumber);
		pointsFile += ".bgeo";

		Utils::writePartio<T, 3>(pointsFile, vertices, vertices->numParticles);
	}

	// helper function
	void inline objToPoly(const std::string &inputFileName, const std::string &outputFileName)
	{
		std::ifstream fin;
		std::ofstream fout;

		fin.open(inputFileName);
		fout.open(outputFileName);

		int vertCtr = 0;
		int faceCtr = 0;
		std::string s;

		std::vector<std::vector<std::string>> vertices;
		std::vector<std::vector<std::string>> faces;

		// READING
        std::string delimiter = "/";

		while(fin >> s) // first one is v, vn or f..
		{
			if (s == "v") // vertex - x, y, z
			{
				std::string x, y, z;
				fin >> x >> y >> z;
                std::vector<std::string> v = {x, y, z};
				vertices.push_back(v);
				vertCtr++;
				// fout << std::to_string(vertCtr) << " " << x << " " << y << " " << z << std::endl;
			}
			else if (s == "vn") // normal - ignore
			{

			}
			else if (s == "f") // face - always a triangle .. v1, v2, v3
			{
				std::string x, y, z;
				fin >> x >> y >> z;

                std::string token = s.substr(0, s.find(delimiter));
                std::vector<std::string> f = {x.substr(0,x.find(delimiter)),
                                              y.substr(0,y.find(delimiter)),
                                              z.substr(0,z.find(delimiter))};
                faces.push_back(f);
				faceCtr++;
				// fout << 1 << std::endl << 3 << " " << x << " " << y << " " << z << std::endl;
			}
			else // comments and stuff
			{

			}
		}
		fin.close();

		// WRITING
        fout << "# OBJ: " << inputFileName << std::endl;
		fout << "# POLY: " << outputFileName << std::endl << std::endl;

        fout << "# vertices" << std::endl;
		fout << std::to_string(vertCtr) << " 3 " << "0 " << "0 " << std::endl << std::endl;
		for(int i=0; i<vertCtr; i++)
		{
			fout << std::to_string(i+1) << " " << vertices[i][0] << " " << vertices[i][1]
									<< " " << vertices[i][2] << std::endl;
		}

		fout << std::endl << "# faces" << std::endl;
		fout << std::to_string(faceCtr) << " 0 " << std::endl << std::endl;
		for(int i=0; i<faceCtr; i++)
		{
			fout << "1 " << std::endl << "3 " << faces[i][0] << " " << faces[i][1]
				 << " " << faces[i][2] << std::endl;
		}

        fout << std::endl << "# holes" << std::endl << "0" << std::endl;
        fout << std::endl << "# regions" << std::endl << "0" << std::endl;

		fout << "# That is all" << std::endl;

		fout.close();
	}

	template<typename T>
	std::vector<T> split(const std::string& line) 
	{
		std::istringstream is(line);
		return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
	}
};
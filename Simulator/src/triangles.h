#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <memory>
#include <stdlib.h>
#include <time.h>

#include "particles.h"
#include "globals.h"
using namespace std;

class Triangles
{
public:
	int numTriangles;
	std::vector<Eigen::Matrix<uint, 3, 1>> triFaceList;
	std::vector<Eigen::Matrix<T, 3, 1>> triNormalList;

	Triangles();
	void addTriangle(int i, int j, int k);

	//create obj file
	void create_objFile(std::string file_name, std::shared_ptr<Particles>& vertices );

	//read tetgen face file
	void tetgen_readLine(std::ifstream &fin);
	void tetgen_readFace(const std::string &inputFileName);

	bool intersect(const Ray &r, const int &triIndex, float *t, std::shared_ptr<Particles>& vertices, Eigen::Matrix<T, 3, 1> *baryCoords) const;
};

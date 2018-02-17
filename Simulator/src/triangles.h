#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

class Triangles
{
public:
	int numTriangles;
	std::vector<Eigen::Matrix<uint,3,1>> triFaceList;

	Triangles();
	void addTriangle(int i, int j, int k);

	//read tetgen face file
	void tetgen_readLine(std::ifstream &fin);
	void tetgen_readFace(const std::string &inputFileName);
};
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <time.h>

class Triangles
{
public:
	int numTriangles;
	std::vector<Eigen::Matrix<uint,3,1>> triFaceList;

public:
	Triangles();
	void addTriangle(int i, int j, int k);
};
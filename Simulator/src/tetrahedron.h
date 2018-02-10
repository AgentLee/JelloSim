#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <time.h>

class Tetrahedron
{
public:
	Tetrahedron();

    Eigen::Matrix<uint, 4, 1> particleIndices;
};
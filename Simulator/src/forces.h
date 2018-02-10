#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <fstream>
using namespace std;

#include "particles.h"
#include "triangles.h"

using T = double;

namespace Forces
{
	Eigen::Matrix<T,3,1> inline getSpringForce( T youngsModulus, T restLength, Eigen::Matrix<T,3,1> pos1, Eigen::Matrix<T,3,1> pos2 )
	{
		Eigen::Matrix<T,3,1> springForce = Eigen::Matrix<T,3,1>::Zero();

		T a = (((pos2-pos1).norm()/restLength) - 1.0f);
		Eigen::Matrix<T,3,1> n12 = (pos1-pos2)/(pos1-pos2).norm();

		springForce = -youngsModulus * a * n12;
		return springForce;
	}

	Eigen::Matrix<T,3,1> inline getSpringDampingForce( T dampingCoeff, Eigen::Matrix<T,3,1> pos1, Eigen::Matrix<T,3,1> pos2, 
		Eigen::Matrix<T,3,1> vel1, Eigen::Matrix<T,3,1> vel2 )
	{
		Eigen::Matrix<T,3,1> dampingForce = Eigen::Matrix<T,3,1>::Zero();

		Eigen::Matrix<T,3,1> n12 = (pos1-pos2)/(pos1-pos2).norm();

		dampingForce = -dampingCoeff * n12 * n12.transpose() * (vel1 - vel2);
		return dampingForce;
	}
};
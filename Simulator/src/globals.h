#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using T = double;

typedef Vector3f Point3f;
typedef Vector4f Point4f;

struct Ray
{
    Eigen::Matrix<T, 3, 1> origin;
    Eigen::Matrix<T, 3, 1> direction;
    Ray(Eigen::Matrix<T, 3, 1> ro, Eigen::Matrix<T, 3, 1> rd): origin(ro), direction(rd) {}
};

struct Intersection
{
	bool hit;
	float t;
	uint triangleIndex;
    Eigen::Matrix<T, 3, 1> point;
    Eigen::Matrix<T, 3, 1> normal;
    Eigen::Matrix<T, 3, 1> BarycentricWeights;// Corresponding to every point that makes the triangle
    Eigen::Matrix<T, 3, 1> penetrationDistance;
    Eigen::Matrix<uint, 3, 1> vertsOfTriangle; //Vertices of the triangle that was intersected
    Eigen::Matrix<T, 3, 1> color;
};

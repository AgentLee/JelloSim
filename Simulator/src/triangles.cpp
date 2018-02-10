#include "triangles.h"

Triangles::Triangles() : numTriangles(0)
{}

void Triangles::addTriangle(int i, int j, int k)
{
	Eigen::Matrix<uint,3,1> tri = {uint(i),uint(j),uint(k)};
	triFaceList.push_back(tri);
	numTriangles++;
}
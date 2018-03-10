#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using T = double;

typedef Vector3f Point3f;
typedef Vector4f Point4f;

struct Ray
{
	Point3f origin;
	Point3f direction;
};

struct Intersection
{
	bool hit;
	float t;
	uint triangleIndex;
	Point3f point;
	Point3f normal;
};

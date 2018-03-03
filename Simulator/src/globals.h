using T = double;
typedef Eigen::Matrix<T, 3, 1> Vector3f;
typedef Eigen::Matrix<T, 4, 1> Vector4f;
typedef Eigen::Matrix<T, 3, 1> Point3f;
typedef Eigen::Matrix<T, 4, 1> Point4f;
typedef Eigen::Matrix<T, 4, 4> Mat4f;

struct Ray
{
	Point3f origin;
	Point3f direction;
};
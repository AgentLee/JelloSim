#include "triangles.h"

Triangles::Triangles() : numTriangles(0)
{}

void Triangles::addTriangle(int i, int j, int k)
{
	Eigen::Matrix<uint, 3, 1> tri = {uint(i),uint(j),uint(k)};
	triFaceList.push_back(tri);
    triNormalList.push_back(Eigen::Matrix<T, 3, 1>::Zero());
	numTriangles++;
}

void Triangles::create_objFile(std::string file_name, std::shared_ptr<Particles>& vertices )
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
    for(int i=0; i<numTriangles; i++)
    {
        objfile << "f  " + std::to_string(triFaceList[i][0] + 1) + "//  "
                         + std::to_string(triFaceList[i][1] + 1) + "//  "
                         + std::to_string(triFaceList[i][2] + 1) + "//  " + "\n";
    }

    objfile << "# end of obj file\n";
    objfile.close();
}

void Triangles::tetgen_readLine(std::ifstream &fin)
{
    float f;
    fin >> f; // first one is id.. //next 3 are the actual indices of the vertices // last one is the boundary marker

    int j, k, l;
    fin >> j >> k >> l;
    addTriangle(j - 1, k - 1, l - 1);

    fin >> f;
}

void Triangles::tetgen_readFace(const std::string &inputFileName)
{
    std::ifstream fin(inputFileName);

    if(fin.is_open())
    {
        int numFaces;
        int boundaryMarker;

        fin >> numFaces >> boundaryMarker;

        for(int i = 0; i < numFaces; ++i)
        {
            tetgen_readLine(fin);
        }

        fin.close();
    }
}

float length(const Eigen::Matrix<T, 3, 1> &v)
{
    return v.norm();
}

Eigen::Matrix<T, 3, 1> cross(const Eigen::Matrix<T, 3, 1> &a, const Eigen::Matrix<T, 3, 1> &b)
{
    return a.cross(b);
}

bool fequal(float a, float b)
{
    return (a < b + 0.001) && (a > b - 0.001);
}

void Triangles::computeNormals(std::shared_ptr<Particles>& vertices)
{
    for(int i=0; i<numTriangles; i++)
    {
        Eigen::Matrix<T, 3, 1> points[3];
        points[0] = vertices->pos[triFaceList[i][0]];
        points[1] = vertices->pos[triFaceList[i][1]];
        points[2] = vertices->pos[triFaceList[i][2]];

        Eigen::Matrix<T, 3, 1> vec1 = points[1] - points[0];
        Eigen::Matrix<T, 3, 1> vec2 = points[2] - points[1];

        triNormalList[i] = -cross(vec1, vec2);
        triNormalList[i].normalize();
    }
}

//The ray in this function is not transformed because it was *already* transformed in Mesh::GetIntersection
bool Triangles::intersect(const Ray &r, const int &triIndex, float *t, std::shared_ptr<Particles>& vertices, Eigen::Matrix<T, 3, 1> *baryCoords) const
{
    Eigen::Matrix<T, 3, 1> planeNormal = triNormalList[triIndex];
    Eigen::Matrix<T, 3, 1> points[3];
    points[0] << vertices->pos[triFaceList[triIndex][0]];
    points[1] << vertices->pos[triFaceList[triIndex][1]];
    points[2] << vertices->pos[triFaceList[triIndex][2]];

//    if(planeNormal.dot(r.direction) < 0)
//    {
//        return false;
//    }

    //1. Ray-plane intersection
    float tnew =  planeNormal.dot(points[0] - r.origin) / planeNormal.dot(r.direction);
    if(tnew < 0.f) return false;

    Eigen::Matrix<T, 3, 1> P = r.origin + tnew * r.direction;
    //2. Barycentric test

    float S = 0.5f * length(cross(points[0] - points[1], points[0] - points[2]));
    float s1 = 0.5f * length(cross(P - points[1], P - points[2]))/S;
    float s2 = 0.5f * length(cross(P - points[2], P - points[0]))/S;
    float s3 = 0.5f * length(cross(P - points[0], P - points[1]))/S;
    float sum = s1 + s2 + s3;

    if(s1 >= 0 && s1 <= 1 && s2 >= 0 && s2 <= 1 && s3 >= 0 && s3 <= 1 && fequal(sum, 1.0f)) {
        *t = tnew;
        (*baryCoords) << s1 , s2, s3;
        return true;
    }

    return false;
}

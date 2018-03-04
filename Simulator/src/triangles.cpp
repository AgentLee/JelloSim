#include "triangles.h"

Triangles::Triangles() : numTriangles(0)
{}

void Triangles::addTriangle(int i, int j, int k)
{
	Eigen::Matrix<uint, 3, 1> tri = {uint(i),uint(j),uint(k)};
	triFaceList.push_back(tri);
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

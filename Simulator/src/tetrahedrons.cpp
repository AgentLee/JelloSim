#include "tetrahedrons.h"
#include <Eigen/Jacobi>

#define oneSixth 0.166666666667

Tetrahedrons::Tetrahedrons()
{
    // TODO
}

void Tetrahedrons::computeRestDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices )
{
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];
    Eigen::Matrix<T, 3, 3> restDef;

    for(uint i=0; i<3; i++)
    {
        //x1-x4  x2-x4  x3-x4
        //y1-y4  y2-y4  y3-y4
        //z1-z4  z2-z4  z3-z4
        restDef(i,0) = vertices->pos[vertexIndices[0]][i] - vertices->pos[vertexIndices[3]][i];
        restDef(i,1) = vertices->pos[vertexIndices[1]][i] - vertices->pos[vertexIndices[3]][i];
        restDef(i,2) = vertices->pos[vertexIndices[2]][i] - vertices->pos[vertexIndices[3]][i];
    }

    restDeformation[tetraIndex] = restDef;
}

void Tetrahedrons::computeInvRestDeformation( uint tetraIndex )
{
    restInverseDeformation[tetraIndex] = restDeformation[tetraIndex].inverse();
}

void Tetrahedrons::computeUndeformedVolume( uint tetraIndex )
{
    undeformedVolume[tetraIndex] = oneSixth * std::abs( restDeformation[tetraIndex].determinant() );
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeNewDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices )
{
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];
    Eigen::Matrix<T, 3, 3> newDef = Eigen::Matrix<T, 3, 3>::Zero();

    for(uint i=0; i<3; i++)
    {
        /*
            x1-x4  x2-x4  x3-x4
            y1-y4  y2-y4  y3-y4
            z1-z4  z2-z4  z3-z4
        */
        newDef(i,0) = vertices->pos[vertexIndices[0]][i] - vertices->pos[vertexIndices[3]][i];
        newDef(i,1) = vertices->pos[vertexIndices[1]][i] - vertices->pos[vertexIndices[3]][i];
        newDef(i,2) = vertices->pos[vertexIndices[2]][i] - vertices->pos[vertexIndices[3]][i];
    }

    return newDef;
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeF( uint tetraIndex, Eigen::Matrix<T,3,3>& Ds )
{
    Eigen::Matrix<T, 3, 3> F = Ds*restDeformation[tetraIndex];
    return F;
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeP( uint tetraIndex, const Eigen::Matrix<T,3,3>& F )
{
    float mu = 0.5;
    float lamda = 0.5;

    Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(F, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix<T, 3, 3> U = svd.matrixU();
    Eigen::Matrix<T, 3, 3> V = svd.matrixV();
    if (U.determinant() < 0.f) {
        U(2) *= -1;
    }
    if (V.determinant() < 0.f) {
        V(2) *= -1;
    }

    Eigen::Matrix<T,3,3> R = U * V.transpose();

    // USE BELOW LINES OF CODE TO USE LINEAR COROTATED METHOD
    // Eigen::Matrix<T,3,3> I = Eigen::Matrix<T, 3, 3>::Identity();
    // Eigen::Matrix<T,3,3> trTerm = R.transpose() * F - I;
    // P = 2.f * mu * (F - R) + lamda * trTerm.trace() * R;

    Eigen::Matrix<T,3,3> J; // jacobian matrix // TODO - Compute J. Eigen doesn't do it??
    double j = J.determinant();
    Eigen::Matrix<T,3,3> JFInvTr  = Eigen::Matrix<T, 3, 3>::Zero();
    if (F.determinant() != 0) {
        JFInvTr = j * F.inverse().transpose(); // TODO - Fix the singular matrix case.
    }

    Eigen::Matrix<T, 3, 3> P = Eigen::Matrix<T, 3, 3>::Zero();
    P = mu * (F - R) + lamda * (j - 1.f) * JFInvTr;
    return P;
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeH( uint tetraIndex, Eigen::Matrix<T,3,3>& P, Eigen::Matrix<T,3,3>& Ds )
{
    Eigen::Matrix<T, 3, 3> H = -(undeformedVolume[tetraIndex] * P * Ds.transpose());
    return H;
}

void Tetrahedrons::addForces( uint tetraIndex, std::shared_ptr<Particles> vertices, Eigen::Matrix<T,3,3>& H )
{
    // Add forces to particles that make up the tetrahedron (f += h)
    // f4 += -(h1 + h2 + h3)
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];
    vertices->force[vertexIndices[0]] += H.col(0);
    vertices->force[vertexIndices[1]] += H.col(1);
    vertices->force[vertexIndices[2]] += H.col(2);
    vertices->force[vertexIndices[3]] += -(H.col(0) + H.col(1) + H.col(2));
}

/*
 *  .node FILE FORMAT:
 *      http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html
 *
 *  First line:
 *              <# of tetrahedra> <nodes per tet. (4 or 10)> <region attribute (0 or 1)>
 *
 *  Remaining lines list # of tetrahedra:
 *              <tetrahedron #> <node> <node> ... <node> [attribute]
 *
 */

// helper function
void Tetrahedrons::tetgen_readLine(std::ifstream &fin, int nodesPerTet)
{
    float f;
    fin >> f; // first one is id..

    for(int i = 0; i < nodesPerTet; ++i)
    {
        fin >> particleIndices[f](i, 0);
    }
}

void Tetrahedrons::tetgen_readEleFile(const std::string &inputFileName)
{
    // TODO
    std::ifstream fin(inputFileName);

    if(fin.is_open())
    {
        int nodesPerTet;
        int region;

        fin >> numTetra >> nodesPerTet >> region;

        particleIndices.resize(numTetra);

        for(int i = 0; i < numTetra; ++i)
        {
            tetgen_readLine(fin, nodesPerTet);
        }

        fin.close();
    }
}

#include "tetrahedrons.h"
#include <Eigen/Jacobi>

#define oneSixth 0.166666666667

Tetrahedrons::Tetrahedrons()
{}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeNewDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices )
{
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];
    Eigen::Matrix<T, 3, 3> newDef = Eigen::Matrix<T, 3, 3>::Zero();

    for(uint i=0; i < 3; i++)
    {
        /*
            x1-x4  x2-x4  x3-x4
            y1-y4  y2-y4  y3-y4
            z1-z4  z2-z4  z3-z4
        */
        newDef(i, 0) = vertices->pos[vertexIndices(0)][i] - vertices->pos[vertexIndices(3)][i];
        newDef(i, 1) = vertices->pos[vertexIndices(1)][i] - vertices->pos[vertexIndices(3)][i];
        newDef(i, 2) = vertices->pos[vertexIndices(2)][i] - vertices->pos[vertexIndices(3)][i];
    }

    return newDef;
}

void Tetrahedrons::computeRestDeformation( uint tetraIndex, std::shared_ptr<Particles> vertices )
{
    restDeformation[tetraIndex] = this->computeNewDeformation( tetraIndex, vertices );
}

void Tetrahedrons::computeInvRestDeformation( uint tetraIndex )
{
    restInverseDeformation[tetraIndex] = restDeformation[tetraIndex].inverse();
}

void Tetrahedrons::computeUndeformedVolume( uint tetraIndex )
{
    undeformedVolume[tetraIndex] = oneSixth * std::abs( restDeformation[tetraIndex].determinant() );
}

void Tetrahedrons::computeUndefVol_into_restInvDefTranspose( uint tetraIndex )
{
    undefVol_into_restInvDefTranspose[tetraIndex] = undeformedVolume[tetraIndex] * restInverseDeformation[tetraIndex].transpose();
}

void Tetrahedrons::addMass( uint tetraIndex, float density, std::shared_ptr<Particles> vertices )
{
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];

    for(uint i=0; i<4; i++)
    {
        vertices->mass[vertexIndices(i)] += 0.25f * density * undeformedVolume[tetraIndex];
    }

}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeF( uint tetraIndex, Eigen::Matrix<T, 3, 3>& Ds )
{
    Eigen::Matrix<T, 3, 3> F = Ds * restInverseDeformation[tetraIndex];

    // HACK to prevent FPN error from building up 
    for(int j = 0; j < 3; ++j) {
        for(int i = 0; i < 3; ++i) {
            if(std::abs(F(i,j)) < 1e-10) {
                F(i, j) = 0;
            } 
        }
    }

    return F;
}

Eigen::Matrix<T, 3, 3> jInvTrMat(const Eigen::Matrix<T,3,3>& mat)
{
    float A, B, C, D, E, F, G, H, I;
    float a, b, c, d, e, f, g, h, i;
    a = mat(0,0), b = mat(0,1), c = mat(0,2),
    d = mat(1,0), e = mat(1,1), f = mat(1,2),
    g = mat(2,0), h = mat(2,1), i = mat(2,2);
    A = e*i-f*h;
    B = f*g-d*i;
    C = d*h-e*g;
    D = c*h-b*i;
    E = a*i-c*g;
    F = b*g-a*h;
    G = b*f-c*e;
    H = c*d-a*f;
    I = a*e-b*d;

    Eigen::Matrix<T, 3, 3> retMat;
    retMat(0,0) = A, retMat(0,1) = B, retMat(0,2) = C,
    retMat(1,0) = D, retMat(1,1) = E, retMat(1,2) = F,
    retMat(2,0) = G, retMat(2,1) = H, retMat(2,2) = I;

    return retMat;
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeP( uint tetraIndex, const Eigen::Matrix<T,3,3>& F, int frame, bool &collided )
{
    // mu = k / (2 * (1 + poisson ratio))
    // lambda = (k * poisson ratio) / ((1 + poisson ratio) (1 - 2 * poisson ratio))
    float youngsMod = 500000.0;
    float poisson = 0.3;   // Make sure this is always less than 0.5 otherwise values go to infinity
    float mu = youngsMod / (2 * (1 + poisson));
    float lamda = (youngsMod * poisson) / ((1 + poisson) * (1 - 2 * poisson));

    Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(F, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix<T, 3, 3> U = svd.matrixU();
    Eigen::Matrix<T, 3, 3> V = svd.matrixV();
    if (U.determinant() < 0.f) {
        U.col(2) *= -1;
    }
    if (V.determinant() < 0.f) {
        V.col(2) *= -1;
    }

    Eigen::Matrix<T, 3, 3> R = U * V.transpose();

    float j = F.determinant();
    Eigen::Matrix<T, 3, 3> JFInvTr  = Eigen::Matrix<T, 3, 3>::Zero();
    JFInvTr = jInvTrMat(F);

    Eigen::Matrix<T, 3, 3> P = Eigen::Matrix<T, 3, 3>::Zero();
    P = 2.0f * mu * (F - R) + lamda * (j - 1.f) * JFInvTr;
    return P;
}

Eigen::Matrix<T, 3, 3> Tetrahedrons::computeH( uint tetraIndex, Eigen::Matrix<T, 3, 3>& P )
{
    Eigen::Matrix<T, 3, 3> H = -( P * undefVol_into_restInvDefTranspose[tetraIndex] );
    return H;
}

void Tetrahedrons::addForces( uint tetraIndex, std::shared_ptr<Particles> vertices, Eigen::Matrix<T, 3, 3>& H )
{
    // Add forces to particles that make up the tetrahedron (f += h)
    // f4 += -(h1 + h2 + h3)
    Eigen::Matrix<uint, 4, 1> vertexIndices = particleIndices[tetraIndex];
    vertices->force[vertexIndices(0)] += H.col(0);
    vertices->force[vertexIndices(1)] += H.col(1);
    vertices->force[vertexIndices(2)] += H.col(2);
    vertices->force[vertexIndices(3)] += -(H.col(0) + H.col(1) + H.col(2));
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
        float num;
        fin >> num;
        particleIndices[f-1](i) = num - 1;
    }
}

void Tetrahedrons::tetgen_readEleFile(const std::string &inputFileName)
{
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

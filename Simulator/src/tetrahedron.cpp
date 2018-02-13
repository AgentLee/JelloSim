#include "tetrahedron.h"

Tetrahedron::Tetrahedron()
{
    // TODO
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
void readLine(std::ifstream &fin, int nodesPerTet)
{
    float f;
    fin >> f; // first one is id..

    for(int i = 0; i < nodesPerTet; ++i)
    {
        fin >> particleIndices[f](i, 0);
    }
}

void Tetrahedron::readEle(const std::string &inputFileName)
{
    // TODO
    std::ifstream fin(inputFileName);

    if(fin.is_open())
    {
        int numTet;
        int nodesPerTet;
        int region;

        fin >> numTet >> nodesPerTet >> region;

        for(int i = 0; i < numTet; ++i)
        {
            readLine(fin, nodesPerTet);
        }

        fin.close();
    }
}

#include "mesh.h"

Mesh::Mesh( float _gridCellSize, float density, float poissonsRatio, float youngsModulus ) 
        : gridCellSize(_gridCellSize)
{
    vertices  = std::make_shared<Particles>(0, 0.0f);
    triangles = std::make_shared<Triangles>();
    tetras    = std::make_shared<Tetrahedrons>();
}

Mesh::Mesh( std::shared_ptr<Particles> _vertices, std::shared_ptr<Triangles> _triangles, std::shared_ptr<Tetrahedrons> _tetras, 
            float _gridCellSize )
        : vertices(_vertices), triangles(_triangles), tetras(_tetras), gridCellSize(_gridCellSize)
{
    //calcBounds
    //create grid
}
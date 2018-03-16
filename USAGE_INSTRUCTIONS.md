# Jello Simulator


## Converting OBJs to Poly Files
Tetgen does not take OBJs as input for tetrahedralization, so we created OBJ to poly converter to load standard models which are available as OBJs. These poly files can then be used with tetgen to generate data that can be used with our code for simulation.

To convert an OBJ to poly, go to `main.cpp` and uncomment `#define CONVERT_OBJ_TO_POLY`

And specify input file name and output file name as follows
```
const std::string objToPolyNames[] = { "../Assets/objs_polys/INPUT.obj", "../Assets/objs_polys/OUTPUT.poly" };
```


## Generating BGEOs for Multiple Files
Our code can simulate multiple simulations simultaneously. 

For simulating multiple meshes, multiple inputs are required. This is achieved as follows

In `main.cpp`
```
// .node files have the vertex data
const std::string nodeFileNames[] = { "../Assets/Meshes/cube_poly_0.001/cube.1.node",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.node" };
// .face files have the data about the triangles of the hull of the meshes
const std::string faceFileNames[] = { "../Assets/Meshes/cube_poly_0.001/cube.1.face",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.face" };
// .ele files have the data of all the tetrahedra in the meshes
const std::string eleFileNames[]  = { "../Assets/Meshes/cube_poly_0.001/cube.1.ele",
									  "../Assets/Meshes/cube_poly_0.5/cube.1.ele"};
// .obj files are generated for visualization with houdini
const std::string objFileNames[]  = { "../Assets/OBJs/FirstCube.obj",
									  "../Assets/OBJs/SecondCube.obj" };
// .bgeo files are used to visualize per point data for every frame in houdini
const std::string bgeoFileNames[] = { "../Assets/BGEOs/jelloCube1Frame",
									  "../Assets/BGEOs/jelloCube2Frame" };
```

Then in `void createScene()` initialize all the meshes with initial translation and other simulation parameters by creating a new shared pointer for Mesh class
```
// mesh 1
{
    std::shared_ptr<Mesh> mesh0 = std::make_shared<Mesh>( nodeFileNames[0], faceFileNames[0], eleFileNames[0], objFileNames[0], gridCellSize, density, poissonsRatio, youngsModulus );
    Vector3f translation = Vector3f(-1.0f, 1.0f, 0.0f);
    mesh0->translateMesh(translation);
    MeshList.push_back(mesh0);
}

// mesh 2
{
    std::shared_ptr<Mesh> mesh1 = std::make_shared<Mesh>( nodeFileNames[1], faceFileNames[1], eleFileNames[1], objFileNames[1], gridCellSize, density, poissonsRatio, youngsModulus );
    Vector3f translation = Vector3f(1.0f, 1.0f, 0.0f);
    mesh1->translateMesh(translation);
    MeshList.push_back(mesh1);
}
```


## Controllable Parameters
`gridCellSize, density, poissonsRatio, and youngsModulus` control the properties of semi-solids in FEM simulation. You can easily change these using globals in `main.cpp`. For a simulation with multiple meshes with different properties, simply set them in `void createScene()` instead of using the globals.

While the inter-object collisions don't work yet, in `sim.cpp`, to toggle interobject collision related settings, change
```
#define SET_POSITIONS 1
#define SET_VELOCITIES 1
#define INTER_OBJECT_COLLISIONS 1
```




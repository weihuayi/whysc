#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/QuadMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QuadMesh;
typedef QuadMesh::Cell Cell;
typedef QuadMesh::Edge Edge;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    QuadMesh mesh;
    MF::one_quad_mesh(mesh);
    mesh.print();

    Writer writer;
    writer.set_points(mesh);
    writer.set_cells(mesh);
    writer.write("test.vtu");

    return 0;
}

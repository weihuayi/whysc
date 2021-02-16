#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;
typedef WHYSC::Mesh::VTKMeshWriter<TriMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    TriMesh mesh;
    MF::square_triangle_mesh(mesh);

    mesh.print();
    Writer writer(&mesh);
    writer.set_points();
    writer.set_cells();
    writer.write("test.vtu");

    return 0;
}

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
typedef WHYSC::Mesh::VTKMeshWriter<HexMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    HexMesh mesh;

    MF::cube_hexahedron_mesh(mesh);
    
    std::vector<double> q;
    mesh.cell_quality(q);

    Writer writer(&mesh);
    writer.set_points();
    writer.set_cells();
    writer.write("test.vtu");

    return 0;
}

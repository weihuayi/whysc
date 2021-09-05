#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
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
    mesh.uniform_refine(5);

    auto NC = mesh.number_of_cells();
    std::vector<int> zdata(NC);

    for(int i = 0; i < NC; i++)
    {
      zdata[i] = mesh.cell_barycenter(i)[2]*10000;
    }
    Writer writer(&mesh);
    writer.set_points();
    writer.set_cells();
    writer.set_cell_data(zdata, 1, "z");
    writer.write("hex_test.vtu");
    return 0;
}

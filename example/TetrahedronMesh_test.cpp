#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    TetMesh mesh;

    //MF::one_tetrahedron_mesh(mesh, "equ");
    //MF::one_tetrahedron_mesh(mesh, "iso");
    MF::cube_tetrahedron_mesh(mesh);
    
    std::vector<double> q;
    mesh.cell_quality(q);
    auto it = std::min_element(q.begin(), q.end());
    std::cout << *it << std::endl;

    for(int i = 0; i < 5; i++)
    {
        double max;
        double min;
        mesh.uniform_refine();
        mesh.cell_quality(q);
        mesh.cell_dihedral_angle(max, min);
        auto it = std::min_element(q.begin(), q.end());
        std::cout << *it << std::endl;
        std::cout << max << ", " << min << std::endl;
    }

    Writer writer;
    writer.set_points(mesh);
    writer.set_cells(mesh);
    writer.write("test.vtu");

    return 0;
}

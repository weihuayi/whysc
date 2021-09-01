#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/IntervalMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::IntervalMesh<GK, Node, Vector> IMesh;
typedef IMesh::Cell Cell;
typedef IMesh::Edge Edge;
typedef WHYSC::Mesh::VTKMeshWriter Writer;

int main(int argc, char **argv)
{
    std::shared_ptr<IMesh> mesh = std::make_shared<IMesh>();
    auto & node = mesh->nodes();
    auto & cell = mesh->cells();

    node.resize(3);
    node[0] = Node(0, 0, 0);
    node[1] = Node(1, 1, 1);
    node[2] = Node(1, 0, 2);
    cell.resize(2);
    cell[0] = Cell({0, 1});
    cell[1] = Cell({1, 2});
    mesh->uniform_refine(4);

    mesh->print();
    Writer writer;
    writer.set_points(*mesh);
    writer.set_cells(*mesh);
    writer.write("test_interval.vtu");
    return 0;
}

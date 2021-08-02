#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"

using namespace WHYSC;

typedef Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;

int main(int argc, char **argv)
{
    TriMesh mesh;

    mesh.insert(Node{0.0, 0.0});
    mesh.insert(Node{1.0, 0.0});
    mesh.insert(Node{1.0, 1.0});
    mesh.insert(Node{0.0, 1.0});

    mesh.insert(Cell{1, 2, 0});
    mesh.insert(Cell{3, 0, 2});
    mesh.init_top(); // 初始化网格拓扑

    std::cout << "Number of nodes:" << mesh.number_of_nodes();
    std::cout << "Number of cells:" << mesh.number_of_cells();
    std::cout << "Number of edges:" << mesh.number_of_edges();
    return 0;
}

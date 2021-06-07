#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"

using namespace WHYSC;

typedef Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;

/*
void test_vector_construction()
{
    std::vector<double> nodes = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
    std::vector<int> cells = {1, 2, 0, 3, 0, 2};

    TriMesh mesh(nodes, cells);

    ASSERT_EQUAL(4, mesh.number_of_nodes());
    ASSERT_EQUAL(2, mesh.number_of_cells());
    ASSERT_EQUAL(5, mesh.number_of_edges());
}
*/

void test_insert_construction()
{
    TriMesh mesh();

    mesh.insert(Node{0.0, 0.0});
    mesh.insert(Node{1.0, 0.0});
    mesh.insert(Node{1.0, 1.0});
    mesh.insert(Node{0.0, 1.0});

    mesh.insert(Cell{1, 2, 0});
    mesh.insert(Cell{3, 0, 2});
    mesh.init_top();
    
    ASSERT_EQUAL(4, mesh.number_of_nodes());
    ASSERT_EQUAL(2, mesh.number_of_cells());
    ASSERT_EQUAL(5, mesh.number_of_edges());
}

int main(int argc, char **argv)
{
  test_insert_construction();
}

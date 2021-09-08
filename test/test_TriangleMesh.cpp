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
    TriMesh mesh;

    mesh.insert(Node{0.0, 0.0});
    mesh.insert(Node{1.0, 0.0});
    mesh.insert(Node{1.0, 1.0});
    mesh.insert(Node{0.0, 1.0});
    mesh.insert(Cell{1, 2, 0});
    mesh.insert(Cell{3, 0, 2});
    mesh.init_top();
    
    //测试节点, 边, 单元数量
    ASSERT_EQUAL(4, mesh.number_of_nodes());
    ASSERT_EQUAL(2, mesh.number_of_cells());
    ASSERT_EQUAL(5, mesh.number_of_edges());

    //测试 edge
    std::vector<std::array<int, 2> > edge_test(5);
    edge_test[0] = std::array<int, 2>({2, 0});
    edge_test[1] = std::array<int, 2>({0, 1});
    edge_test[2] = std::array<int, 2>({1, 2});
    edge_test[3] = std::array<int, 2>({2, 3});
    edge_test[4] = std::array<int, 2>({3, 0});
    for(int i = 0; i < 5; i++)
    {
      ASSERT_EQUAL(edge_test[i][0], mesh.edge(i)[0]);
      ASSERT_EQUAL(edge_test[i][1], mesh.edge(i)[1]);
    }

    //测试面积
    ASSERT_THROW(mesh.cell_measure(0)-0.5<1e-10);
    ASSERT_THROW(mesh.cell_measure(1)-0.5<1e-10);

    //测试 edge2cell
    std::vector<std::array<int, 4> > edge2cell_test(5);
    edge2cell_test[0] = std::array<int, 4>({0, 1, 0, 0});
    edge2cell_test[1] = std::array<int, 4>({0, 0, 1, 1});
    edge2cell_test[2] = std::array<int, 4>({0, 0, 2, 2});
    edge2cell_test[3] = std::array<int, 4>({1, 1, 1, 1});
    edge2cell_test[4] = std::array<int, 4>({1, 1, 2, 2});
    for(int i = 0; i < 5; i++)
    {
      ASSERT_EQUAL(edge2cell_test[i][0], mesh.edge_to_cell(i)[0]);
      ASSERT_EQUAL(edge2cell_test[i][1], mesh.edge_to_cell(i)[1]);
      ASSERT_EQUAL(edge2cell_test[i][2], mesh.edge_to_cell(i)[2]);
      ASSERT_EQUAL(edge2cell_test[i][3], mesh.edge_to_cell(i)[3]);
    }
}

int main(int argc, char **argv)
{
  test_insert_construction();
}

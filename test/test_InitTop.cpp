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

void test_insert_construction()
{
    TriMesh mesh;
    int N = 100;
    auto & nodes = mesh.nodes();
    auto & cells = mesh.cells();

    cells.resize(N);
    nodes.resize(3*N);
    for(int i = 0; i < N; i++)
    {
        nodes[3*i] = Node{0.0, 0.0};
        nodes[3*i+1] = Node{1.0, 0.0};
        nodes[3*i+2] = Node{0.0, 1.0};
        cells[i] = Cell({3*i, 3*i+1, 3*i+2});
    }
    auto start = clock();
    mesh.init_top0();
    auto end = clock();
    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;
}

int main(int argc, char **argv)
{
  test_insert_construction();
}

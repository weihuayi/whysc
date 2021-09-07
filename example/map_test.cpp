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
typedef TetMesh::Cell Cell;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    int N = std::stoi(argv[1]);
    TetMesh mesh;
    auto & nodes = mesh.nodes();
    auto & cells = mesh.cells();

    cells.resize(N);
    nodes.resize(4*N);
    for(int i = 0; i < N; i++)
    {
        nodes[4*i] = Node{0.0, 0.0, 0.0};
        nodes[4*i+1] = Node{1.0, 0.0, 0.0};
        nodes[4*i+2] = Node{0.0, 1.0, 0.0};
        nodes[4*i+3] = Node{0.0, 0.0, 1.0}; 
        cells[i] = Cell({4*i, 4*i+1, 4*i+2, 4*i+3});
    }
    auto start = clock();
    mesh.init_top();
    auto end = clock();
    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;

    TetMesh mesh0;
    auto & node0s = mesh0.nodes();
    auto & cell0s = mesh0.cells();

    cell0s.resize(N);
    node0s.resize(4*N);
    for(int i = 0; i < N; i++)
    {
        node0s[4*i] = Node{0.0, 0.0, 0.0};
        node0s[4*i+1] = Node{1.0, 0.0, 0.0};
        node0s[4*i+2] = Node{0.0, 1.0, 0.0};
        node0s[4*i+3] = Node{0.0, 0.0, 1.0}; 
        cell0s[i] = Cell({4*i, 4*i+1, 4*i+2, 4*i+3});
    }
    start = clock();
    mesh0.init_top0();
    end = clock();
    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;
    return 0;
}

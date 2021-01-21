#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;

int main(int argc, char **argv)
{
    TriMesh mesh;
    // 插入四个点
    mesh.insert(Node{0.0, 0.0});
    mesh.insert(Node{1.0, 0.0});
    mesh.insert(Node{1.0, 1.0});
    mesh.insert(Node{0.0, 1.0});

    mesh.insert(Cell{1, 2, 0});
    mesh.insert(Cell{3, 0, 2});
    
    mesh.number_of_holes() = 1; // 上面有网格有一个外部无限大的洞
    mesh.number_of_genus() = 0; // 亏格为 0

    mesh.construct_top();
    mesh.print();
    mesh.uniform_refine();

    auto NC = mesh.number_of_cells();
    for(int i = 0; i < NC; i++)
    {
        std::cout << i << ": " << mesh.cell_measure(i) << std::endl;
    }

    std::vector<double> eh;
    mesh.edge_measure(eh);
    std::cout << eh.size() << std::endl;
    for(int i = 0; i < eh.size(); i++)
    {
        std::cout << i << ": " <<  eh[i]<< std::endl;
    }



    return 0;
}

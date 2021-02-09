#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef TetMesh::Cell Cell;

int main(int argc, char **argv)
{
    TetMesh mesh;

    mesh.insert(Node{0.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 1.0, 1.0});
    mesh.insert(Node{0.0, 1.0, 1.0});

    mesh.insert(Cell{0, 1, 2, 6});
    mesh.insert(Cell{0, 5, 1, 6});
    mesh.insert(Cell{0, 4, 5, 6});
    mesh.insert(Cell{0, 7, 4, 6});
    mesh.insert(Cell{0, 3, 7, 6});
    mesh.insert(Cell{0, 2, 3, 6});
    
    mesh.init_top();
    mesh.print();

    mesh.uniform_refine();
    mesh.print();

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

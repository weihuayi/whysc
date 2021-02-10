#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/VTKMeshWriter.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef TetMesh::Cell Cell;
typedef WHYSC::Mesh::VTKMeshWriter<TetMesh> Writer;

int main(int argc, char **argv)
{
    TetMesh mesh;
/*
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
*/

    mesh.insert(Node{0.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 0.0, 0.0});
    mesh.insert(Node{0.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 0.0, 1.0});

    mesh.insert(Cell{0, 1, 2, 3});
    
    mesh.init_top();
    
    std::vector<double> q;
    mesh.cell_quality(q);
    auto it = std::min_element(q.begin(), q.end());
    std::cout << *it << std::endl;

    for(int n = 0; n < 6; n++)
    {
        mesh.uniform_refine();
        mesh.cell_quality(q);
        auto it = std::min_element(q.begin(), q.end());
        std::cout << *it << std::endl;
    }

    Writer writer(&mesh);
    writer.write();

    return 0;
}

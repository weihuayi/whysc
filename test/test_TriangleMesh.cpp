
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::Mesh::TriangleMesh<GK> Mesh;
typedef Mesh::

int main(int argc, char **argv)
{
    std::vector<double> nodes = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
    std::vector<int> cells = {1, 2, 0, 3, 0, 2};

    Mesh tri(nodes, cells);
    tri.print();

    tri.uniform_refine(nodes, cells);

    Mesh newTri(nodes, cells);
    return 0;
}

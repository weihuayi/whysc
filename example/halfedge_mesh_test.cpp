#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/HalfEdgeMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::HalfEdgeMesh<GK, Node, Vector> HMesh;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char **argv)
{
    TriMesh tmesh;
    MF::one_triangle_mesh(tmesh);
    tmesh.uniform_refine(1);

    HMesh hmesh;
    hmesh.from_mesh(tmesh);
    hmesh.print();
    return 0;
}

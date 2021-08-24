#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/HexJacobiQuality.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
typedef WHYSC::Mesh::HexJacobiQuality<HexMesh> MeshQuality;
typedef WHYSC::Mesh::MeshFactory MF;
typedef HexMesh::Vector Vector;
typedef HexMesh::Toplogy Toplogy;

int main()
{
  auto mesh = std::make_shared<HexMesh>();
  MF::cube_hexahedron_mesh(*mesh);
  //mesh->uniform_refine();
  auto & node = mesh->nodes();
  node[0] =Node(-1, -1, -1);
  node[0] =Node(0.2, 0.2, 0.2);
  for(auto n : node)
  {
    std::cout<< "n = " << n <<std::endl;
  }

  MeshQuality tmq(mesh);
  int NC = mesh->number_of_cells();
  for(int i = 0; i < NC; i++)
  {
    auto v = tmq.quality(i);
    std::cout<< "v = " << v <<std::endl;
  }
  auto qqq = tmq.gradient(0, 0);
  std::cout<<  " gradient = " << qqq <<std::endl;
  return 0;
}

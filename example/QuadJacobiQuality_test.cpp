#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/QuadMesh.h"
#include "mesh/QuadJacobiQuality.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
typedef WHYSC::Mesh::QuadJacobiQuality<QMesh> MeshQuality;
typedef WHYSC::Mesh::MeshFactory MF;
typedef QMesh::Vector Vector;
typedef QMesh::Toplogy Toplogy;

int main()
{
  auto mesh = std::make_shared<QMesh>();
  MF::one_quad_mesh(*mesh);
  //mesh->uniform_refine();
  auto & node = mesh->nodes();
  //node[3][2] = 10;
  node[0] =Node(0.2, 0.2);
  for(auto n : node)
  {
    std::cout<< "n = " << n[0] << " " << n[1] <<std::endl;
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

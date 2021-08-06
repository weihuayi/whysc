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
  for(auto n : node)
  {
    std::cout<< "n = " << n[0] << " " << n[1] <<std::endl;
  }

  MeshQuality tmq(mesh);
  Vector v;
  double w;
  tmq.nabla_quality(0, 3, v, w);
  std::cout<< "lll - "<<tmq.quality_of_cell(0) <<std::endl;
  std::cout<< "v = " << v[0] << " " << v[1] << " " << v[2] <<std::endl;
  std::cout<< "w = " << w <<std::endl;

  /*
  Toplogy n2c;
  mesh->node_to_cell(n2c);
  for(auto q0 : n2c.local_indices())
    std::cout<< "q = " << q0 <<std::endl;
  for(auto q0 : n2c.locations())
    std::cout<< "loc = " << q0 <<std::endl;
    */
  return 0;
}

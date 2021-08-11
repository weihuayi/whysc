#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/TriRadiusRatioQuality.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::TriRadiusRatioQuality<TriMesh> TriMeshQuality;
typedef WHYSC::Mesh::MeshFactory MF;
typedef TriMesh::Vector Vector;
typedef TriMesh::Toplogy Toplogy;

int main()
{
  auto mesh = std::make_shared<TriMesh>();
  MF::one_triangle_mesh(*mesh, "equ");
  //mesh->uniform_refine();
  auto & node = mesh->nodes();
  //node[2][1] = 10;
  node[2][0] = 0;
  node[2][1] = 1;
  node[0][0] = 0.1;
  node[0][1] = 0.1;
  node[0]=Node(0.11844, 0.46333)+(Node(0, 0)-Node(0.11844, 0.46333));
  node[1]=Node(0.08451, 0.37679)+(Node(0, 0)-Node(0.11844, 0.46333));
  node[2]=Node(0.14119, 0.40032)+(Node(0, 0)-Node(0.11844, 0.46333));
  node[0]=Node(0, 0);
  node[1]=Node(-6.13022832e-02, -6.98745059e-02);
  node[2]=Node( 4.55805753e-18, -6.69912129e-02);
  //node[0]=Node(0, 0);
  //node[1]=Node(-1, -1.13);
  //node[2]=Node(0, -1.09);
  for(auto n : node)
  {
    std::cout<< "n = " << n[0] << " " << n[1] <<std::endl;
  }

  TriMeshQuality tmq(mesh);
  auto v = tmq.gradient(0, 2);
  auto q = tmq.quality(0);
  std::cout<< "cell_q = " << q <<std::endl;
  std::cout<< "v = " << v[0] << " " << v[1]<<std::endl;

  node[0]=Node(0, 0);
  node[1]=Node(-1, -1.13);
  node[2]=Node(0, -1.09);
  for(auto n : node)
  {
    std::cout<< "n = " << n[0] << " " << n[1] <<std::endl;
  }

  v = tmq.gradient(0, 2);
  q = tmq.quality(0);
  std::cout<< "cell_q = " << q <<std::endl;
  std::cout<< "v = " << v[0] << " " << v[1]<<std::endl;

  Toplogy n2c;
  mesh->node_to_cell(n2c);
  for(auto q0 : n2c.local_indices())
    std::cout<< "q = " << q0 <<std::endl;
  for(auto q0 : n2c.locations())
    std::cout<< "loc = " << q0 <<std::endl;
  return 0;
}

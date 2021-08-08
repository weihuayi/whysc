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
  auto & cell = mesh->cells();
  //node[3][2] = 10;
  for(auto n : node)
  {
    std::cout<< "n = " << n[0] << " " << n[1] <<std::endl;
  }

  std::vector<const Node *> nodeArray;
  for(int i = 0; i < 4; i++)
  {
    nodeArray.push_back(&(node[cell[0][i]]));
    std::cout<< "lx: " << *nodeArray.back() <<std::endl;
  }

  MeshQuality tmq;
  auto v = tmq.nabla(nodeArray);
  auto qqq = tmq.quality(nodeArray);
  std::cout<< "cell_q = " << qqq <<std::endl;
  std::cout<< "v = " << v[0] << " " << v[1] <<std::endl;

  auto & idx = mesh->m_num[3];
  for(int j = 0; j < 4; j++)
  {
    nodeArray[j] = &(mesh->node(idx[j]));
  }
  v = tmq.nabla(nodeArray);
  std::cout<< "v = " << v[0] << " " << v[1] <<std::endl;
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

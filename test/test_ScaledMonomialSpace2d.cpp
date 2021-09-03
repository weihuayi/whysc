#include <math.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"
#include "mesh/QuadMesh.h"
#include "mesh/TriangleMesh.h"
#include "functionspace/ScaledMonomialSpace2d.h"

typedef WHYSC::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;


void test_basis_on_quad_mesh(int p=1)
{
  typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
  typedef QMesh::Cell Cell;
  typedef QMesh::Edge Edge;
  typedef WHYSC::FunctionSpace::ScaledMonomialSpace2d<QMesh, AK> Space;

  auto mesh = std::make_shared<QMesh>();
  auto & nodes = mesh->nodes();
  auto & cells = mesh->cells();

  nodes.resize(6);
  cells.resize(2);
  nodes[0] = Node({0.0, 0.0});
  nodes[1] = Node({2.0, 0.0});
  nodes[2] = Node({2.0, 2.0});
  nodes[3] = Node({0.0, 2.0});
  nodes[4] = Node({4.0, 0.0});
  nodes[5] = Node({4.0, 2.0});
  cells[0] = Cell({0, 1, 2, 3});
  cells[1] = Cell({1, 4, 5, 2});
  mesh->init_top();

  std::vector<Node> point(2);
  point[0] = Node({0, 0.5});
  point[1] = Node({1, 0.5});

  auto space = std::make_shared<Space>(mesh, p);

  Matrix val(2, (p+1)*(p+2)/2);
  space->basic(point, val);
  std::cout<< val <<std::endl;
}

int main(int argc, char **argv)
{
  test_basis_on_quad_mesh();
}

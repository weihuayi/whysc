#include <math.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"

#include "mesh/QuadMesh.h"

#include "functionspace/ScaledMonomialSpace2d.h"

typedef WHYSC::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
typedef QMesh::Cell Cell;
typedef QMesh::Edge Edge;
typedef WHYSC::FunctionSpace::ScaledMonomialSpace2d<QMesh, AK> Space;

int main(int argc, char **argv)
{
  auto mesh = std::make_shared<QMesh>();
  auto & node = mesh->nodes();
  auto & cell = mesh->cells();

  node.resize(6);
  cell.resize(2);
  node[0] = Node({0.0, 0.0});
  node[1] = Node({2.0, 0.0});
  node[2] = Node({2.0, 2.0});
  node[3] = Node({0.0, 2.0});
  node[4] = Node({4.0, 0.0});
  node[5] = Node({4.0, 2.0});
  cell[0] = Cell({0, 1, 2, 3});
  cell[1] = Cell({1, 4, 5, 2});
  mesh->init_top();

  std::vector<Node> point(2);
  point[0] = Node({0, 0.5});
  point[1] = Node({1, 0.5});

  int p = std::stoi(argv[1]);
  auto space = std::make_shared<Space>(mesh, p);

  Matrix val;
  space->basic(point, val);
  std::cout<< val <<std::endl;

  Matrix grad_xval;
  Matrix grad_yval;
  space->grad_basic(point, grad_xval, grad_yval);
  std::cout<< grad_xval <<std::endl;
  std::cout<< grad_yval <<std::endl;

  Matrix lval;
  space->laplace_basic(point, lval);
  std::cout<< lval <<std::endl;
}













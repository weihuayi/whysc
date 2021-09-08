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
#include "mesh/TetrahedronMesh.h"
#include "mesh/HexahedronMesh.h"
#include "functionspace/ScaledMonomialSpace2d.h"

typedef WHYSC::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;

double monomial(int idx, const Node & point)
{
  int p = 0;
  for(int i = 0; i < 10; i++)
  {
    if(idx > i)
    {
      idx -= i+1;
      p += 1;
    }
    else
      break;
  }
  return std::pow(point[0], p-idx)*std::pow(point[1], idx);
}

Vector grad_monomial(int idx, const Node & point)
{
  int p = 0;
  for(int i = 0; i < 10; i++)
  {
    if(idx > i)
    {
      idx -= i+1;
      p += 1;
    }
    else
      break;
  }
  auto s0 = (p-idx)*std::pow(point[0], p-idx-1)*std::pow(point[1], idx);
  auto s1 = idx*std::pow(point[0], p-idx)*std::pow(point[1], idx-1);
  return Vector({s0, s1});
}

double laplace_monomial(int idx, const Node & point)
{
  int p = 0;
  for(int i = 0; i < 10; i++)
  {
    if(idx > i)
    {
      idx -= i+1;
      p += 1;
    }
    else
      break;
  }
  auto s0 = (p-idx-1)*(p-idx)*std::pow(point[0], p-idx-2)*std::pow(point[1], idx);
  auto s1 = (idx-1)*idx*std::pow(point[0], p-idx)*std::pow(point[1], idx-2);
  return s0+s1;
}

void test_basis_on_quad_mesh(int p=1)
{
  typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
  typedef QMesh::Cell Cell;
  typedef QMesh::Edge Edge;
  typedef WHYSC::FunctionSpace::ScaledMonomialSpace2d<QMesh> Space;

  auto mesh = std::make_shared<QMesh>();
  auto & nodes = mesh->nodes();
  auto & cells = mesh->cells();

  nodes.resize(4);
  cells.resize(1);
  nodes[0] = Node({0.0, 0.0});
  nodes[1] = Node({2.0, 0.0});
  nodes[2] = Node({2.0, 2.0});
  nodes[3] = Node({0.0, 2.0});
  cells[0] = Cell({0, 1, 2, 3});
  mesh->init_top();

  auto space = std::make_shared<Space>(mesh, p);

  Node point = Node({0, 0.5});
  Node center;
  mesh->cell_barycenter(0, center);
  auto h = mesh->cell_size(0);

  std::vector<double> val;
  space->basis(0, point, val);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = monomial(i, (point-center)/h);
    ASSERT_THROW(s-val[i]<1e-10);
  }

  std::vector<Vector> gval;
  space->grad_basis(0, point, gval);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = grad_monomial(i, (point-center)/h);
    ASSERT_THROW(s[0]-gval[i][0]<1e-10);
    ASSERT_THROW(s[1]-gval[i][1]<1e-10);
  }

  std::vector<double> lval;
  space->laplace_basis(0, point, lval);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = laplace_monomial(i, (point-center)/h);
    ASSERT_THROW(s-lval[i]<1e-10);
  }
}

void test_basis_on_tri_mesh(int p=1)
{
  typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> Mesh;
  typedef Mesh::Cell Cell;
  typedef Mesh::Edge Edge;
  typedef WHYSC::FunctionSpace::ScaledMonomialSpace2d<Mesh> Space;

  auto mesh = std::make_shared<Mesh>();
  auto & nodes = mesh->nodes();
  auto & cells = mesh->cells();

  nodes.resize(4);
  cells.resize(1);
  nodes[0] = Node({0.0, 0.0});
  nodes[1] = Node({2.0, 0.0});
  nodes[2] = Node({2.0, 2.0});
  cells[0] = Cell({0, 1, 2});
  mesh->init_top();

  auto space = std::make_shared<Space>(mesh, p);

  Node point = Node({0, 0.5});
  Node center;
  mesh->cell_barycenter(0, center);
  auto h = mesh->cell_size(0);

  std::vector<double> val;
  space->basis(0, point, val);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = monomial(i, (point-center)/h);
    ASSERT_THROW(s-val[i]<1e-10);
  }

  std::vector<Vector> gval;
  space->grad_basis(0, point, gval);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = grad_monomial(i, (point-center)/h);
    ASSERT_THROW(s[0]-gval[i][0]<1e-10);
    ASSERT_THROW(s[1]-gval[i][1]<1e-10);
  }

  std::vector<double> lval;
  space->laplace_basis(0, point, lval);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    auto s = laplace_monomial(i, (point-center)/h);
    ASSERT_THROW(s-lval[i]<1e-10);
  }
}

int main(int argc, char **argv)
{
  test_basis_on_quad_mesh(9);
  test_basis_on_tri_mesh(9);
}

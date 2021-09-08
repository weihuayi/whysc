#ifndef QuadPositiveJacobiQuality_h
#define QuadPositiveJacobiQuality_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>
#include <algorithm>
#include <limits.h>

#include "CellQualityBase.h"

namespace WHYSC {
namespace Mesh {

/*
 *  
 * Note
 * ----
 * 四面体外接球与内接球半径之比
 *       mu = R/r/3
 * 
 * 给定四面体网格, 可以计算这个网格任意一个 patch 的质量.
 *
*/

template<typename Mesh>
class QuadPositiveJacobiQuality: public CellQualityBase<Mesh>
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

public:
  using CellQualityBase<Mesh>::CellQualityBase;

  /*
   * 网格第 i 个单元的质量
   */
  double quality(int i);

  /*
   * 网格第 c 个单元的质量关于这个单元的第 i 个点的梯度
   */
  Vector gradient(int c, int i);

};

template<typename Mesh>
inline double QuadPositiveJacobiQuality<Mesh>::quality(int i)
{
  auto mesh = this->get_mesh();
  auto & cell = mesh->cell(i);
  auto & node = mesh->nodes();

  bool flag = false;
  double q = 0.0;
  for(int j = 0; j < 4; j++)
  {
    auto current = cell[j];

    auto v1 = node[cell[(j+1)%4]] - node[current];
    auto v2 = node[cell[(j+3)%4]] - node[current];

    auto L1 = v1.squared_length();
    auto L2 = v2.squared_length();

    auto J = cross(v1, v2);
    q += (L1+L2)/(2.0*J);
    if(J<0)
    {
      flag = true;
      break;
    }
  }
  if(flag)
  {
    return LONG_MAX;
  }
  else
  {
    return q/4.0;
  }
}

template<typename Mesh>
inline typename Mesh::Vector QuadPositiveJacobiQuality<Mesh>::gradient(int c, int i)
  //第 c 个单元的质量对单元的第 i 个点求梯度
{
  auto mesh = this->get_mesh();
  auto & node = mesh->nodes();
  auto & cell = mesh->cell(c);

  //四个顶点, n1, n2, n3 是逆时针
  auto p0 = cell[i];
  auto p1 = cell[(i+1)%4];
  auto p2 = cell[(i+2)%4];
  auto p3 = cell[(i+3)%4];

  auto v01 = node[p1] - node[p0];
  auto v12 = node[p2] - node[p1];
  auto v23 = node[p3] - node[p2];
  auto v30 = node[p0] - node[p3];

  //三条边的长度
  auto L1 = v01.squared_length();
  auto L2 = v12.squared_length();
  auto L3 = v23.squared_length();
  auto L4 = v30.squared_length();

  //计算 nabla_mu0
  auto J = cross(v30, v01);
  auto d = L1+L4;
  auto q = d/(2.0*J);

  auto nabla_J = Vector(v01[1]+v30[1], -v30[0]-v01[0]);
  auto nabla_d = 2.0*(v01-v30);
  auto nabla_q0 = (-d*nabla_J + J*nabla_d)/(2.0*J*J);

  //计算 nabla_mu1
  J = cross(v01, v12);
  d = L2+L1;
  q = d/(2.0*J);

  nabla_J = Vector(-v12[1], v12[0]);
  nabla_d = -2.0*v01;
  auto nabla_q1 = (-d*nabla_J + J*nabla_d)/(2.0*J*J);

  //计算 nabla_mu3
  J = cross(v23, v30);
  d = L4+L3;
  q = d/(2.0*J);

  nabla_J = Vector(-v23[1], v23[0]);
  nabla_d = 2.0*v30;
  auto nabla_q3 = (-d*nabla_J + J*nabla_d)/(2.0*J*J);

  auto v = nabla_q0 + nabla_q1 + nabla_q3;
  return v;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadPositiveJacobiQuality_h

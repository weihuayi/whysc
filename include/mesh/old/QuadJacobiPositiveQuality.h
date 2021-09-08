#ifndef QuadJacobiPositiveQuality_h
#define QuadJacobiPositiveQuality_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>

namespace WHYSC {
namespace Mesh {

/*
 *  
 * Note
 * ----
 * 
 * 给定四面体网格, 可以计算这个网格任意一个 patch 的质量.
 */
template<typename QMesh>
class QuadJacobiPositiveQuality
{
public:
  typedef typename QMesh::Node Node;
  typedef typename QMesh::Vector Vector;

public:
  QuadJacobiPositiveQuality(std::shared_ptr<QMesh> mesh)
  {
    m_mesh = mesh;
  }

  // 计算第 i 个单元的质量
  double quality_of_cell(int i)
  {
    auto & cell = m_mesh->cells();
    auto & node = m_mesh->nodes();
    double q = 0;

    for(int j = 0; j < 4; j++)
    {
      auto current = cell[i][j];

      auto v1 = node[cell[i][(j+1)%4]] - node[current];
      auto v2 = node[cell[i][(j+3)%4]] - node[current];

      auto L1 = v1.squared_length();
      auto L2 = v2.squared_length();

      auto J = cross(v1, v2);
      q += (L1+L2)/(2.0*J);
    }
    return q;
  }

  template<typename Patch>
  double patch_quality(Patch & patch)
  {
    double q = 0;
    int N = patch.number_of_adj_entities();
    for(int i = 0; i < N; i++)
    {
      q += quality_of_cell(patch.adj_entity(i));
    }
    return q;
  }

  void nabla_quality(int c, int i, Vector &v, double &w)//第 c 个单元的质量对单元的第 i 个点求梯度
  {
    auto & node = m_mesh->nodes();
    auto & cell = m_mesh->cell(c);

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

    v = nabla_q0 + nabla_q1 + nabla_q3;
    v *= -1.0;
    w = 1.0;
  }

private:
  std::shared_ptr<QMesh> m_mesh;
};


} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadJacobiPositiveQuality_h

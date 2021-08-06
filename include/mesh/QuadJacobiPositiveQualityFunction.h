#ifndef QuadJacobiPositiveQualityFunction_h
#define QuadJacobiPositiveQualityFunction_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>
#include "BaseObjectFunction.h"

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
class QuadJacobiPositiveQualityFunction: public BaseObjectionFunction<QMesh> 
{
public:
  typedef typename QMesh::Node Node;
  typedef typename QMesh::Vector Vector;
  typedef typename QMesh::Toplogy Toplogy;

public:
  QuadJacobiPositiveQualityFunction(std::shared_ptr<QMesh> mesh): 
     BaseObjectionFunction<QMesh>(mesh)
  {
    mesh->node_to_cell(m_n2c);
  }
  // 计算第 i 个单元的质量
  double quality_of_cell(int i)
  {
    auto mesh = this->get_mesh();
    auto & cell = mesh->cells();
    auto & node = mesh->nodes();
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

  double get_value(int i)
  {
    double q = 0;
    auto patch = m_n2c.adj_entities_with_local(i);

    int N = patch.number_of_adj_entities();
    for(int i = 0; i < N; i++)
    {
      q += quality_of_cell(patch.adj_entity(i));
    }
    return q;
  }

  
  void get_grad_value(int i, Vector & v)
  {
    auto patch = m_n2c.adj_entities_with_local(i);
    int NP = patch.number_of_adj_entities();

    Vector v0 = {0};
    for(int k = 0; k < NP; k++)
    {
      Vector tmpVector;
      nabla(patch.adj_entity(k), patch.adj_local_index(k), tmpVector);
      v0 = v0 + tmpVector;
    }
    double w = min_len(i);
    v = w*v0/std::sqrt(v0.squared_length()); //std::sqrt(v0.squared_length());
  }

  void nabla(int c, int i, Vector &v)//第 c 个单元的质量对单元的第 i 个点求梯度
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

    v = nabla_q0 + nabla_q1 + nabla_q3;
    v *= -1.0;
  }

  double min_len(int i)
  {
    auto mesh = this->get_mesh();
    auto patch = m_n2c.adj_entities_with_local(i);
    auto & cell = mesh->cells();
    auto & node = mesh->nodes();
    int NP = patch.number_of_adj_entities();
    double L = 100000.0;

    for(int i = 0; i < NP; i++)
    {
      auto & cell = mesh->cell(patch.adj_entity(i));
      auto j = patch.adj_local_index(i);
      auto current = cell[j];

      auto v1 = node[cell[(j+1)%4]] - node[current];
      double L1 = std::sqrt(v1.squared_length());
      if(L>L1)
        L=L1;
    }
    return L;
  }


private:
  Toplogy m_n2c;
};


} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadJacobiPositiveQualityFunction_h

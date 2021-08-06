#ifndef TetRadiusRatioQualityFunction_h
#define TetRadiusRatioQualityFunction_h

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
 * 四面体外接球与内接球半径之比
 *       mu = R/r/3
 * 
 * 给定四面体网格, 可以计算这个网格任意一个 patch 的质量.
 *
*/

template<typename TMesh>
class TetRadiusRatioQualityFunction: public BaseObjectionFunction<TMesh> 
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;
  typedef typename TMesh::Toplogy Toplogy;

public:
  // 计算第 i 个单元的质量
  TetRadiusRatioQualityFunction(std::shared_ptr<TMesh> mesh): 
     BaseObjectionFunction<TMesh>(mesh)
  {
    mesh->node_to_cell(m_n2c);
  }

  double quality_of_cell(int i)
  {
    auto mesh = this->get_mesh();
    auto & cell = mesh->cells();
    auto & node = mesh->nodes();

    auto v1 = node[cell[i][1]] - node[cell[i][0]];
    auto v2 = node[cell[i][2]] - node[cell[i][0]];
    auto v3 = node[cell[i][3]] - node[cell[i][0]];

    auto L1 = v1.squared_length();
    auto L2 = v2.squared_length();
    auto L3 = v3.squared_length();

    auto S12 = std::sqrt(cross(v1, v2).squared_length())/2;
    auto S23 = std::sqrt(cross(v2, v3).squared_length())/2;
    auto S31 = std::sqrt(cross(v3, v1).squared_length())/2;
    auto S123 = std::sqrt(cross(v3-v1, v2-v1).squared_length())/2;

    auto S = S12+S23+S31+S123;
    auto d = L1*cross(v2, v3) + L2*cross(v3, v1) + L3*cross(v1, v2);
    auto Ld = std::sqrt(d.squared_length());
    auto V = dot(cross(v1, v2), v3)/6;
    auto r = 3*V/S;
    auto R = Ld/V/12;
    return R/r/3;
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

    double w = 0;
    Vector v0 = {0};
    for(int k = 0; k < NP; k++)
    {
      Vector tmpVector;
      double tmpW;
      nabla(patch.adj_entity(k), patch.adj_local_index(k), tmpVector, tmpW);
      v0 = v0 + tmpVector;
      w = w + tmpW;
    }
    v = v0/w; //std::sqrt(v0.squared_length());
  }

  void nabla(int c, int i, Vector &v, double &w)//第 c 个单元的质量对单元的第 i 个点求梯度
  {
    auto & index = TMesh::m_localface;
    auto mesh = this->get_mesh();
    auto & node = mesh->nodes();
    auto & cell = mesh->cell(c);

    //四个顶点, n1, n2, n3 是逆时针
    auto n0 = cell[i];
    auto n1 = cell[index[i][0]];
    auto n2 = cell[index[i][1]];
    auto n3 = cell[index[i][2]];

    //除去底面外的三条边的向量
    auto v1 = node[n1] - node[n0];
    auto v2 = node[n2] - node[n0];
    auto v3 = node[n3] - node[n0];

    //三条边的长度
    auto L1 = v1.squared_length();
    auto L2 = v2.squared_length();
    auto L3 = v3.squared_length();

    //三个向量之间的互相叉乘
    auto c12 = cross(v1, v2);
    auto c23 = cross(v2, v3);
    auto c31 = cross(v3, v1);

    auto S12 = std::sqrt(c12.squared_length())/2;
    auto S23 = std::sqrt(c23.squared_length())/2;
    auto S31 = std::sqrt(c31.squared_length())/2;
    auto S123 = std::sqrt(cross(v3-v1, v2-v1).squared_length())/2;

    auto S = S12+S23+S31+S123;
    auto d = L1*c23 + L2*c31 + L3*c12;
    auto Ld = std::sqrt(d.squared_length());
    auto V = dot(c12, v3)/6;
    auto r = 3*V/S;
    auto R = Ld/V/12;
    auto Q = R/r/3;

    auto nabla_Ld = ((2.0*(dot(d, c23)*v1 + dot(d, c31)*v2 + dot(d, c12)*v3))
      + cross(d, (L3-L2)*v1 + (L1-L3)*v2 + (L2-L1)*v3))/Ld;
    auto nabla_s = (cross(c12, v1-v2)/S12 + cross(c23, v2-v3)/S23 + cross(c31, v3-v1)/S31)/4.0;
    auto nabla_V = (c12 + c23 + c31)/6.0;

    auto w_d = 2*dot(d, c12+c23+c31)/d.squared_length();
    auto w_s = ((v1-v2).squared_length()/S12 + (v2-v3).squared_length()/S23 + 
      (v3-v1).squared_length()/S31)/S/4.0;

    v = Q*(nabla_Ld/Ld + nabla_s/S - 2.0*nabla_V/V);
    w = w_d+w_s;
  }

private:
  Toplogy m_n2c;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of TetRadiusRatioQualityFunction_h

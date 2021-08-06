#ifndef TriRadiusRatioQualityFunction_h
#define TriRadiusRatioQualityFunction_h

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
class TriRadiusRatioQualityFunction: public BaseObjectionFunction<TMesh> 
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;
  typedef typename TMesh::Toplogy Toplogy;

public:
  // 计算第 i 个单元的质量
  TriRadiusRatioQualityFunction(std::shared_ptr<TMesh> mesh): 
     BaseObjectionFunction<TMesh>(mesh)
  {
    mesh->node_to_cell(m_n2c);
  }

  double quality_of_cell(int i);
  /*
   * 计算第 i 个单元的质量
   */

  double get_value(int i);
  /*
   * 得到第 i 的节点周围单元的质量和
   */

  void get_grad_value(int i, Vector & v);
  /*
   * 得到得到总质量的第 i 的节点的偏导
   */

  void nabla(int c, int i, Vector &v);
  /*
   * 得到第 c 个单元的质量对单元的第 i 个点求梯度
   */

  double min_len(int i);
  /*
   * 得到连接第 i 个点的最短边长度
   */

private:
  Toplogy m_n2c;
};


template<typename TMesh>
inline double TriRadiusRatioQualityFunction<TMesh>::quality_of_cell(int i)
{
  auto mesh = this->get_mesh();
  auto & cell = mesh->cells();
  auto & node = mesh->nodes();

  auto v1 = node[cell[i][1]] - node[cell[i][0]];
  auto v2 = node[cell[i][2]] - node[cell[i][0]];
  auto v12 = node[cell[i][2]] - node[cell[i][1]];

  auto L1 = std::sqrt(v1.squared_length());
  auto L2 = std::sqrt(v2.squared_length());
  auto L12 = std::sqrt(v12.squared_length());

  auto A = cross(v1, v2)/2;
  auto R = L1*L2*L12/A/4;
  auto r = 2*A/(L1+L2+L12);
  return R/r/2;
}

template<typename TMesh>
inline double TriRadiusRatioQualityFunction<TMesh>::get_value(int i)
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

template<typename TMesh>
inline void TriRadiusRatioQualityFunction<TMesh>::get_grad_value(int i, Vector & v)
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

template<typename TMesh>
inline void TriRadiusRatioQualityFunction<TMesh>::nabla(int c, int i, 
    Vector &v)//第 c 个单元的质量对单元的第 i 个点求梯度
{
  auto mesh = this->get_mesh();
  auto & node = mesh->nodes();
  auto & cell = mesh->cell(c);

  //四个顶点, n1, n2, n3 是逆时针
  auto x0 = cell[i];
  auto x1 = cell[(i+1)%3];
  auto x2 = cell[(i+2)%3];

  //除去底面外的三条边的向量
  auto v1 = node[x1] - node[x0];
  auto v2 = node[x2] - node[x0];
  auto v12 = node[x2] - node[x1];

  //三条边的长度
  auto l1 = v1.squared_length();
  auto l2 = v2.squared_length();
  auto l12 = v12.squared_length();

  auto A = cross(v1, v2)/2.0;
  auto R = l1*l2*l12/A/4.0;
  auto r = 2.0*A/(l1+l2+l12);

  auto nabla_l1 = -v1/l1;
  auto nabla_l2 = -v2/l2;
  auto nabla_A = Vector(-v12[1], v12[0])/2.0;
  auto nabla_R = l12*(A*(l2*nabla_l1 + l1*nabla_l2) - l1*l2*nabla_A)/A/A/4.0;
  auto nabla_r = 2.0*((l1+l2+l12)*nabla_A - A*(nabla_l1+nabla_l2))/(l1+l2+l12)/(l1+l2+l12);
  auto nabla_q = (r*nabla_R - R*nabla_r)/r/r/2.0;
  v = nabla_q;
}

template<typename TMesh>
inline double TriRadiusRatioQualityFunction<TMesh>::min_len(int i)
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

    auto v1 = node[cell[(j+1)%3]] - node[current];
    double L1 = std::sqrt(v1.squared_length());
    if(L>L1)
      L=L1;
  }
  return L;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of TriRadiusRatioQualityFunction_h

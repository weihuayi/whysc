#ifndef HexJacobiQualityFunction_h
#define HexJacobiQualityFunction_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>
#include <algorithm>

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
class HexJacobiQualityFunction: public BaseObjectionFunction<TMesh> 
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;
  typedef typename TMesh::Toplogy Toplogy;

public:
  // 计算第 i 个单元的质量
  HexJacobiQualityFunction(std::shared_ptr<TMesh> mesh): 
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
inline double HexJacobiQualityFunction<TMesh>::quality_of_cell(int i)
{
  auto mesh = this->get_mesh();
  const auto & cell = mesh->cell(i);
  auto & node = mesh->nodes();
  auto & index = mesh->m_num;

  double mu = 0.0;
  for(auto & idx : index)
  {
    const auto & x0 = node[cell[idx[0]]];
    const auto & x4 = node[cell[idx[4]]];
    const auto & x2 = node[cell[idx[2]]];
    const auto & x1 = node[cell[idx[1]]];

    auto v1 = x4 - x0;
    auto v2 = x2 - x0;
    auto v3 = x1 - x0;
    
    auto l1 = std::sqrt(v1.squared_length());
    auto l2 = std::sqrt(v1.squared_length());
    auto l3 = std::sqrt(v1.squared_length());

    auto J = dot(cross(v1, v2), v3);
    auto d = l1*l1*l1 + l2*l2*l2 + l3*l3*l3;
    auto q = 3*J/d;
    mu += std::pow((q-2), 4);
  }
  return mu;
}

template<typename TMesh>
inline double HexJacobiQualityFunction<TMesh>::get_value(int i)
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
inline void HexJacobiQualityFunction<TMesh>::get_grad_value(int i, Vector & v)
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
inline void HexJacobiQualityFunction<TMesh>::nabla(int c, int i,  Vector &v)
  //第 c 个单元的质量对单元的第 i 个点求梯度
{
  auto mesh = this->get_mesh();
  const auto & node = mesh->nodes();
  const auto & cell = mesh->cell(c);
  const auto & idx = mesh->m_num[i];

  const auto & x0 = node[cell[idx[0]]]; 
  const auto & x1 = node[cell[idx[1]]]; 
  const auto & x2 = node[cell[idx[2]]]; 
  const auto & x3 = node[cell[idx[3]]]; 
  const auto & x4 = node[cell[idx[4]]]; 
  const auto & x5 = node[cell[idx[5]]]; 
  const auto & x6 = node[cell[idx[6]]]; 

  auto v04 = x4-x0;  
  auto v02 = x2-x0;  
  auto v01 = x1-x0;  
  auto v46 = x6-x4;  
  auto v45 = x5-x4;  
  auto v23 = x3-x2;  
  auto v26 = x6-x2;  
  auto v13 = x3-x1;  
  auto v15 = x5-x1;  

  auto l04 = std::sqrt(v04.squared_length());  
  auto l02 = std::sqrt(v02.squared_length());  
  auto l01 = std::sqrt(v01.squared_length());  
  auto l46 = std::sqrt(v46.squared_length());  
  auto l45 = std::sqrt(v45.squared_length());  
  auto l23 = std::sqrt(v23.squared_length());  
  auto l26 = std::sqrt(v26.squared_length());  
  auto l13 = std::sqrt(v13.squared_length());  
  auto l15 = std::sqrt(v15.squared_length());  

  auto J0 = dot(cross(v04, v02), v01);
  auto J4 = dot(cross(v46, v45), v04);
  auto J2 = dot(cross(v23, v26), v02);
  auto J1 = dot(cross(v15, v13), v01);

  auto d0 = std::pow(l04, 3) + std::pow(l01, 3) + std::pow(l02, 3);
  auto d4 = std::pow(l04, 3) + std::pow(l45, 3) + std::pow(l46, 3);
  auto d2 = std::pow(l26, 3) + std::pow(l23, 3) + std::pow(l02, 3);
  auto d1 = std::pow(l13, 3) + std::pow(l01, 3) + std::pow(l15, 3);

  auto q0 = 3*J0/d0;  
  auto q4 = 3*J4/d4;  
  auto q2 = 3*J2/d2;  
  auto q1 = 3*J1/d1;  

  auto nabla_J0 = corss(v04, v01) + cross(v02, v04) + cross(v01, v02);
  auto nabla_J4 = cross(v45, v46);
  auto nabla_J2 = cross(v26, v23);
  auto nabla_J1 = cross(v13, v15);

  auto nabla_d0 = -3*(l04*v04 + l02*v02 + l01*v01);
  auto nabla_d4 = 3*l04*v04;
  auto nabla_d2 = 3*l02*v02;
  auto nabla_d1 = 3*l01*v01;

  auto nabla_q0 = 3*(d0*nabla_J0 - J0*nabla_d0)/d0/d0;
  auto nabla_q4 = 3*(d4*nabla_J4 - J4*nabla_d4)/d4/d4;
  auto nabla_q2 = 3*(d2*nabla_J2 - J2*nabla_d2)/d2/d2;
  auto nabla_q1 = 3*(d1*nabla_J1 - J1*nabla_d1)/d1/d1;

  auto nabla_mu0 = 4*std::pow(q0 - 2, 3)*nabla_q0;
  auto nabla_mu4 = 4*std::pow(q4 - 2, 3)*nabla_q4;
  auto nabla_mu2 = 4*std::pow(q2 - 2, 3)*nabla_q2;
  auto nabla_mu1 = 4*std::pow(q1 - 2, 3)*nabla_q1;
  v = nabla_mu0 + nabla_mu4 + nabla_mu2 + nabla_mu1;
}

template<typename TMesh>
inline double HexJacobiQualityFunction<TMesh>::min_len(int i)
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
    auto & idx = mesh->m_num[j];
    auto x0 = node[cell[idx[0]]];
    auto x4 = node[cell[idx[4]]];
    auto x2 = node[cell[idx[2]]];
    auto x1 = node[cell[idx[1]]];

    double l[3];
    l[0] = std::sqrt((x4-x0).squared_length());
    l[1] = std::sqrt((x2-x0).squared_length());
    l[2] = std::sqrt((x1-x0).squared_length());

    double L1 = *std::min_element(l, l+3);
    if(L>L1)
      L=L1;
  }
  return L;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of HexJacobiQualityFunction_h

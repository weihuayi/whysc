#ifndef TetRadiusRatioQuality_h
#define TetRadiusRatioQuality_h

#include <vector>
#include <memory>
#include <math.h>

#include "CellQualityBase.h"

namespace WHYSC {
namespace Mesh {

/*
 *  
 * Note
 * ----
 * 四面体外接球与内接球半径之比
 *       mu = R/r/3
 */

template<typename TMesh>
class TetRadiusRatioQuality: public CellQualityBase<TMesh> 
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;

public:
  using CellQualityBase<TMesh>::CellQualityBase;

  double quality(int i);
  /*
   * 网格第 i 个单元的质量
   */


  Vector gradient(int c, int i);
  /*
   * 网格第 c 个单元的质量关于这个单元的第 i 个点的梯度
   */

};

template<typename TMesh>
inline double TetRadiusRatioQuality<TMesh>::quality(int i)
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

template<typename TMesh>
inline typename TMesh::Vector TetRadiusRatioQuality<TMesh>::gradient(int c, int i)
{
  auto mesh = this->get_mesh();
  auto & index = TMesh::m_localface;
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

  auto v = -Q*(nabla_Ld/Ld + nabla_s/S - 2.0*nabla_V/V);
  //w = w_d+w_s;
  return v;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of TetRadiusRatioQuality_h

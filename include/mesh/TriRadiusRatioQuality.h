#ifndef TriRadiusRatioQuality_h
#define TriRadiusRatioQuality_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>

#include "CellQualityBase.h"

namespace WHYSC {
namespace Mesh {

/*
 *  
 * Note
 * ----
 * 三角形外接圆与内接圆半径之比
 *       mu = R/r/2
*/

template<typename TMesh>
class TriRadiusRatioQuality: public CellQualityBase<TMesh>
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
inline double TriRadiusRatioQuality<TMesh>::quality(int i)
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
inline typename TMesh::Vector TriRadiusRatioQuality<TMesh>::gradient(int c, int i)//第 c 个单元的质量对单元的第 i 个点求梯度
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
  auto l1 = std::sqrt(v1.squared_length());
  auto l2 = std::sqrt(v2.squared_length());
  auto l12= std::sqrt(v12.squared_length());

  auto A = cross(v1, v2)/2.0;
  auto R = l1*l2*l12/A/4.0;
  auto r = 2.0*A/(l1+l2+l12);

  auto nabla_l1 = -v1/l1;
  auto nabla_l2 = -v2/l2;
  auto nabla_A = Vector(-v12[1], v12[0])/2.0;
  auto nabla_R = l12*(A*(l2*nabla_l1 + l1*nabla_l2) - l1*l2*nabla_A)/A/A/4.0;
  auto nabla_r = 2.0*((l1+l2+l12)*nabla_A - A*(nabla_l1+nabla_l2))/(l1+l2+l12)/(l1+l2+l12);
  auto nabla_q = (r*nabla_R - R*nabla_r)/r/r/2.0;
  return nabla_q;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of TriRadiusRatioQuality_h

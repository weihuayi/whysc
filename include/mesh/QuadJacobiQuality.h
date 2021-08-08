#ifndef QuadJacobiQuality_h
#define QuadJacobiQuality_h

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
 * 四边形单元质量
 */
template<typename QMesh>
class QuadJacobiQuality: public CellQualityBase<QMesh>
{
public:
  typedef typename QMesh::Node Node;
  typedef typename QMesh::Vector Vector;

public:
  double quality(const std::vector<const Node*> & vertex);
  /*
   * 以 `vertex` 为顶点的单元的质量
   */

  Vector nabla(const std::vector<const Node*> & vertex);
  /*
   * 以 `vertex` 为顶点的单元的质量对第 0 个顶点的梯度
   */

};

template<typename QMesh>
inline double QuadJacobiQuality<QMesh>::quality(const std::vector<const Node*> & vertex)
  {
    double q = 0;
    for(int j = 0; j < 4; j++)
    {
      auto v1 = *(vertex[(j+1)%4]) - *(vertex[j]);
      auto v2 = *(vertex[(j+3)%4]) - *(vertex[j]);

      auto L1 = v1.squared_length();
      auto L2 = v2.squared_length();

      auto J = cross(v1, v2);
      auto mu = 2.0*J/(L1+L2);
      q += std::pow(mu-2.0, 4);
    }
    return q;
  }

template<typename QMesh>
inline typename QMesh::Vector QuadJacobiQuality<QMesh>::nabla(const std::vector<const Node*> & vertex)
  {
    auto v01 = *vertex[1] - *vertex[0];
    auto v12 = *vertex[2] - *vertex[1];
    auto v23 = *vertex[3] - *vertex[2];
    auto v30 = *vertex[0] - *vertex[3];

    //三条边的长度
    auto L1 = v01.squared_length();
    auto L2 = v12.squared_length();
    auto L3 = v23.squared_length();
    auto L4 = v30.squared_length();

    //计算 nabla_mu0
    auto J = cross(v30, v01);
    auto d = L1+L4;
    auto q = 2.0*J/d;

    auto nabla_J = Vector(v01[1]+v30[1], -v30[0]-v01[0]);
    auto nabla_d = 2.0*(v01-v30);
    auto nabla_q = (2.0/(d*d))*(d*nabla_J - J*nabla_d);
    auto nabla_mu0 = 4.0*std::pow(q-2.0, 3)*nabla_q;

    //计算 nabla_mu1
    J = cross(v01, v12);
    d = L2+L1;
    q = 2.0*J/d;

    nabla_J = Vector(-v12[1], v12[0]);
    nabla_d = -2.0*v01;
    nabla_q = (2.0/(d*d))*(d*nabla_J - J*nabla_d);
    auto nabla_mu1 = 4.0*std::pow(q-2.0, 3)*nabla_q;

    //计算 nabla_mu3
    J = cross(v23, v30);
    d = L4+L3;
    q = 2.0*J/d;

    nabla_J = Vector(-v23[1], v23[0]);
    nabla_d = 2.0*v30;
    nabla_q = (2.0/(d*d))*(d*nabla_J - J*nabla_d);
    auto nabla_mu3 = 4.0*std::pow(q-2, 3)*nabla_q;

    return -(nabla_mu0 + nabla_mu1 + nabla_mu3);
  }

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadJacobiQuality_h

#ifndef QuadJacobiQuality_h
#define QuadJacobiQuality_h

#include <vector>
#include <array>
#include <map>
#include <memory>
#include <math.h>
#include <limits.h>

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
  using CellQualityBase<QMesh>::CellQualityBase;

  /*
   * 网格第 i 个单元的质量
   */
  double quality(int i);

  /*
   * 网格第 c 个单元的质量关于这个单元的第 i 个点的梯度
   */
  Vector gradient(int c, int i);

};

template<typename QMesh>
inline double QuadJacobiQuality<QMesh>::quality(int i)
  {
    auto mesh = this->get_mesh();
    auto & cell = mesh->cells();
    auto & node = mesh->nodes();
    double q = 0;

    int flag = 0;
    for(int j = 0; j < 4; j++)
    {
      auto current = cell[i][j];

      auto v1 = node[cell[i][(j+1)%4]] - node[current];
      auto v2 = node[cell[i][(j+3)%4]] - node[current];

      auto L1 = v1.squared_length();
      auto L2 = v2.squared_length();

      auto J = cross(v1, v2);
      auto mu = 2.0*J/(L1+L2);
      q += std::pow(mu-2.0, 4);
      if(J < 0){flag+=1;}
    }

    if(flag >1)
    {
      return LONG_MAX;
    }
    else
    {
      return q/4.0;
    }
  }

template<typename QMesh>
inline typename QMesh::Vector QuadJacobiQuality<QMesh>::gradient(int c, int i)
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

    return nabla_mu0 + nabla_mu1 + nabla_mu3;
  }

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadJacobiQuality_h

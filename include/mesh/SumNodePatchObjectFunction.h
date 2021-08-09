#ifndef SumNodePatchObjectFunction_h
#define SumNodePatchObjectFunction_h

/*
 * 文件: 优化过程中函数对象的基类
*/

#include <math.h>
#include <memory>
#include <vector>

#include "NodePatchObjectFunctionBase.h"

namespace WHYSC{
namespace Mesh{

template<typename Mesh, typename CellQuality>
class SumNodePatchObjectFunction: public NodePatchObjectFunctionBase<Mesh, CellQuality>
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

public:
  SumNodePatchObjectFunction(const std::shared_ptr<Mesh> mesh):
    NodePatchObjectFunctionBase<Mesh, CellQuality>(mesh){}

  double value(Vector); 
  /* 
   * 计算当 m_patch 中心为 node + Vector时, patch 的质量
   */

  Vector gradient(); 
  /*
   * 计算当前 patch 的质量关于 patch 中心节点的梯度.
   */
};


template<typename Mesh, typename CellQuality>
inline double SumNodePatchObjectFunction<Mesh, CellQuality>::value(Vector v)
{
  auto mesh = this->get_mesh();
  auto & patch = this->get_patch();
  auto pid = patch.id();
  auto cell_quality = this->get_cell_quality();

  mesh->node(pid) = mesh->node(pid) + v;

  double q = 0;
  int N = patch.number_of_adj_entities();
  for(int i = 0; i < N; i++)
  {
    q += cell_quality->quality(patch.adj_entity(i)); 
  }
  mesh->node(pid) = this->get_node();
  return q;
}

template<typename Mesh, typename CellQuality>
inline typename Mesh::Vector SumNodePatchObjectFunction<Mesh, CellQuality>::gradient() 
{
  auto mesh = this->get_mesh();
  auto cell_quality = this->get_cell_quality();
  auto & patch = this->get_patch();

  Vector v = {0};

  int NP = patch.number_of_adj_entities();
  for(int i = 0; i < NP; i++)
  {
    v = v + cell_quality->gradient(patch.adj_entity(i), patch.adj_local_index(i)); 
  }
  auto w = this->min_len()/3;
  v = w*v/std::sqrt(v.squared_length());
  return v;
}


};//end of Mesh
};//end of WHYSC

#endif // end of SumNodePatchObjectFunction_h

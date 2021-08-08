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

  double value(Node); 
  /* 
   * 计算当 m_patch 中心为 node 时, patch 的质量
   */

  Vector gradient(); 
  /*
   * 计算当前 patch 的质量关于 patch 中心节点的梯度.
   */
};


template<typename Mesh, typename CellQuality>
inline double SumNodePatchObjectFunction<Mesh, CellQuality>::value(Node node)
{
  auto mesh = this->get_mesh();
  auto & patch = this->get_patch();
  auto cell_quality = this->get_cell_quality();
  int NNC = mesh->number_of_nodes_of_each_cell();

  std::vector<const Node* > nodeArray(NNC);
  nodeArray[0] = &node;

  double q = 0;
  int N = patch.number_of_adj_entities();
  for(int i = 0; i < N; i++)
  {
    auto & cell = mesh->cell(patch.adj_entity(i));
    auto & idx = mesh->m_num[patch.adj_local_index(i)];
    for(int j = 1; j < NNC; j++)
    {
      nodeArray[j] = &(mesh->node(cell[idx[j]]));
    }
    q += cell_quality->quality(nodeArray); 
  }
  return q;
}

template<typename Mesh, typename CellQuality>
inline typename Mesh::Vector SumNodePatchObjectFunction<Mesh, CellQuality>::gradient() 
{
  auto mesh = this->get_mesh();
  auto cell_quality = this->get_cell_quality();
  auto & patch = this->get_patch();
  int NNC = mesh->number_of_nodes_of_each_cell();

  Vector v = {0};
  std::vector<const Node* > nodeArray(NNC);

  int NP = patch.number_of_adj_entities();
  for(int i = 0; i < NP; i++)
  {
    auto & cell = mesh->cell(patch.adj_entity(i));
    auto & idx = mesh->m_num[patch.adj_local_index(i)];
    for(int j = 0; j < NNC; j++)
    {
      nodeArray[j] = &(mesh->node(cell[idx[j]]));
    }
    v = v + cell_quality->nabla(nodeArray); 
  }
  auto w = this->min_len()/3;
  v = w*v/std::sqrt(v.squared_length());
  return v;
}


};//end of Mesh
};//end of WHYSC

#endif // end of SumNodePatchObjectFunction_h

#ifndef NodePatchObjectFunctionBase_h
#define NodePatchObjectFunctionBase_h

/*
 * 文件: 优化过程中函数对象的基类
*/

#include <cmath>
#include <map>
#include <memory>
#include <vector>
#include <math.h>
#include <algorithm>

namespace WHYSC{
namespace Mesh{

template<typename Mesh, typename CellQuality>
class NodePatchObjectFunctionBase
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;
  typedef typename Mesh::Toplogy Toplogy;
  typedef typename Mesh::Toplogy::AdjEntitySetWithLoc Patch;

public:
  NodePatchObjectFunctionBase(std::shared_ptr<Mesh> mesh, 
                              Patch * patch):
    m_mesh(mesh), m_patch(patch)
  {
    m_node = mesh->node(patch->id());
    m_cell_quality = std::make_shared<CellQuality>(mesh);
  }

  /* 
   * 计算当 m_patch 中心为 node + Vector时, patch 的质量
   */
  virtual double value(Node) = 0; 

  /*
   * 计算当前 patch 的质量关于 patch 中心节点的梯度.
   */
  virtual Vector gradient() = 0; 

  virtual Vector direction() = 0;

  Node & get_node() {return m_node;}

  Patch * get_patch() {return m_patch;}

  std::shared_ptr<Mesh> get_mesh() {return m_mesh;}

  std::shared_ptr<CellQuality> get_cell_quality() {return m_cell_quality;}

  double min_len()
  {
    auto mesh = this->get_mesh();
    auto patch = this->get_patch();

    auto & cell = mesh->cells();
    auto & node = mesh->nodes();

    int NP = patch->number_of_adj_entities();
    int NNC = mesh->number_of_nodes_of_each_cell();

    double L = 100000.0;

    for(int i = 0; i < NP; i++)
    {
      auto & cell = mesh->cell(patch->adj_entity(i));
      auto j = patch->adj_local_index(i);
      auto & idx = mesh->m_num[j];
      for(int k = 1; k < NNC; k++)
      {
        auto v = node[cell[idx[k]]] - node[cell[idx[0]]];
        auto L1 = std::sqrt(v.squared_length());
        if(L>L1)
          L=L1;
      }
    }
    return L;
  }

private:
  Node m_node;
  Patch * m_patch;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<CellQuality> m_cell_quality;
};




};//end of Mesh
};//end of WHYSC

#endif // end of NodePatchObjectFunctionBase_h

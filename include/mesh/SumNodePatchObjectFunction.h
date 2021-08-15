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
  typedef typename Mesh::Toplogy Toplogy;
  typedef typename Mesh::Toplogy::AdjEntitySetWithLoc Patch;


public:
  using NodePatchObjectFunctionBase<Mesh, CellQuality>::NodePatchObjectFunctionBase;

  SumNodePatchObjectFunction(std::shared_ptr<Mesh> mesh, Patch * patch):
    m_mesh(mesh), m_patch(patch)
  {
    m_node = mesh->node(patch->id());
    m_cell_quality = std::make_shared<CellQuality>(mesh);
  }

  /* 
   * 计算当 m_patch 中心为 node + Vector时, patch 的质量
   */
  double value(Node); 

  /*
   * 计算当前 patch 的质量关于 patch 中心节点的梯度.
   */
  Vector gradient(); 

  Vector direction();

  Patch * get_patch()
  {
    return m_patch;
  }

private:
  Node m_node;
  Patch * m_patch;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<CellQuality> m_cell_quality;
};

template<typename Mesh, typename CellQuality>
inline double SumNodePatchObjectFunction<Mesh, CellQuality>::value(Node n)
{
  auto pid = m_patch->id();
  m_mesh->node(pid) = n;

  double q = 0;
  int N = m_patch->number_of_adj_entities();
  for(int i = 0; i < N; i++)
  {
    q += m_cell_quality->quality(m_patch->adj_entity(i)); 
  }
  m_mesh->node(pid) = m_node;
  return q;
}

template<typename Mesh, typename CellQuality>
inline typename Mesh::Vector SumNodePatchObjectFunction<Mesh, CellQuality>::gradient() 
{
  Vector v = {0};

  int NP = m_patch->number_of_adj_entities();
  for(int i = 0; i < NP; i++)
  {
    v = v + m_cell_quality->gradient(m_patch->adj_entity(i), m_patch->adj_local_index(i)); 
  }
  return v;
}

template<typename Mesh, typename CellQuality>
inline typename Mesh::Vector SumNodePatchObjectFunction<Mesh, CellQuality>::direction() 
{
  return -gradient();
}

};//end of Mesh
};//end of WHYSC

#endif // end of SumNodePatchObjectFunction_h

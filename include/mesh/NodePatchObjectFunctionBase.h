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

public:
  NodePatchObjectFunctionBase() {} 

  /* 
   * 计算当 m_patch 中心为 node + Vector时, patch 的质量
   */
  virtual double value(Node) = 0; 

  /*
   * 计算当前 patch 的质量关于 patch 中心节点的梯度.
   */
  virtual Vector gradient() = 0; 

  virtual Vector direction() = 0;
};




};//end of Mesh
};//end of WHYSC

#endif // end of NodePatchObjectFunctionBase_h

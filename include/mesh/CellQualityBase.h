#ifndef CellQualityBase_h
#define CellQualityBase_h

#include <vector>
namespace WHYSC {
namespace Mesh {

template<typename TMesh>
class CellQualityBase
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;

public:
  // 计算第 i 个单元的质量
  CellQualityBase(){} 

  virtual double quality(const std::vector<const Node*> & ) = 0;

  virtual Vector nabla(const std::vector<const Node*> & ) = 0;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of CellQualityBase_h

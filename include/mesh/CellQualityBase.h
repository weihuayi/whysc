#ifndef CellQualityBase_h
#define CellQualityBase_h

#include <vector>
#include <memory>

namespace WHYSC {
namespace Mesh {

template<typename TMesh>
class CellQualityBase
{
public:
  typedef typename TMesh::Node Node;
  typedef typename TMesh::Vector Vector;

public:
  CellQualityBase(std::shared_ptr<TMesh> mesh): m_mesh(mesh){} 

  virtual double quality(int i) = 0;
  /*
   * 网格第 i 个单元的质量
   */

  virtual Vector gradient(int c, int i) = 0;
  /*
   * 网格第 c 个单元的质量关于这个单元的第 i 个点的梯度
   */

  void quality_of_mesh(std::vector<double> & q) // 所有单元的质量
  {
    int NC = m_mesh->number_of_cells();
    q.resize(NC);
    for(int i = 0; i < NC; i++)
    {
      q[i] = quality(i);
    }
  }

  std::shared_ptr<TMesh> get_mesh()
  {
    return m_mesh;
  }

private:
  std::shared_ptr<TMesh> m_mesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of CellQualityBase_h

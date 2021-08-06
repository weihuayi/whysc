#ifndef BaseObjectionFunction_h
#define BaseObjectionFunction_h

/*
 * 文件: 优化过程中函数对象的基类
*/

#include <map>
#include <memory>
#include <vector>

namespace WHYSC{
namespace Mesh{

template<typename M>
class BaseObjectionFunction
{
public:
  typedef typename M::Vector Vector;
public:
  BaseObjectionFunction(const std::shared_ptr<M> mesh): m_mesh(mesh){}

  virtual double get_value(int){return 0;} 

  virtual void get_grad_value(int, Vector){}; 

  std::shared_ptr<M> get_mesh()
  {
    return m_mesh;
  }

  double id()
  {
    return 0.1;
  }

private:
  std::shared_ptr<M> m_mesh;
};




};//end of Mesh
};//end of WHYSC

#endif // end of BaseObjectionFunction_h

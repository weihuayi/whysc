#ifndef QuadratureBase_h
#define QuadratureBase_h

#include <vector>
#include <memory>

namespace WHYSC {
namespace Quadrature {

class QuadratureBase
{
public:
  QuadratureBase(){} 

  virtual int number_of_quadrature_points(); 

  virtual void get_quadrature_point_and_weight();
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of QuadratureBase_h

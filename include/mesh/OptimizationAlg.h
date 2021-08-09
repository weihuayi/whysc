#ifndef OptimizationAlg_h
#define OptimizationAlg_h

/*
 * 文件: 各种优化程序的库文件.
*/

#include <math.h>
#include <functional>

namespace WHYSC{
namespace Mesh{

class OptimizationAlg
{

public:
  static double line_search(std::function<double(double)> F);
};

inline double OptimizationAlg::line_search(std::function<double(double)> F)
{
  double a = 0;
  double b = 1;
  double c = a + (1 - 0.618)*(b-a);
  double d = a + 0.618*(b-a);

  double Qc = F(c);
  double Qd = F(d);

  while(abs(Qc-Qd)>0.0001)// 0.618法
  {
    if(Qc > Qd)
    {
      a = c;
      c = d;
      d = a + 0.618*(b-a);
      Qc = Qd;
      Qd = F(d);
    }
    else
    {
      b = d;
      d = c;
      c = a + (1 - 0.618)*(b-a);
      Qd = Qc;
      Qc = F(c);
    }
  }
  return c;
}


};//end of Mesh
};//end of WHYSC

#endif // end of OptimizationAlg_h

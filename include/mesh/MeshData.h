#ifndef MeshData_h
#define MeshData_h

#include <vector>

namespace WHYSC {
namespace Mesh {

struct NodeData
{
  NodeData(){}
  void resize(int n)
  {
    gdof.resize(n);
    gtag.resize(n);
  }
  std::vector<int> gdof;
  std::vector<int> gtag;
};

}
}

#endif // end of MeshData_h

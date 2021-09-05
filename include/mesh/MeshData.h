#ifndef MeshData_h
#define MeshData_h

#include <vector>

namespace WHYSC {
namespace Mesh {

struct NodeData
{
  NodeData(){}
  NodeData(std::vector<int> idx, NodeData data)
  {
    reinit(idx, data);
  }

  void reinit(std::vector<int> idx, NodeData data)
  {
    int N = idx.size();
    resize(N);
    for(int i = 0; i < N; i++)
    {
      gdof[i] = data.gdof[idx[i]];
      gtag[i] = data.gtag[idx[i]];
    }
  }

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

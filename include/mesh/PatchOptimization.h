#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg;

template<typename PMesh>
class PatchOptimization
{
public:
  typedef typename PMesh::Node Node;
  typedef typename PMesh::Toplogy Toplogy;

public:
  PatchOptimization(std::shared_ptr<PMesh> mesh)
  {
    m_mesh = mesh;
    mesh->node_to_cell(m_top);
  }

  void optimization(int i, Node & node)
  {
    int GD = m_mesh->geo_dimension();
    for(int n = 0; n < GD; n++)
      node[n] = 0;

    int NP = patch_size(i);
    for(int k = 0; k < NP; k++)
    {
      Node tmpNode;
      m_mesh->cell_barycenter(patch(i, k), tmpNode);
      for(int j = 0; j < GD; j++)
      {
        node[j] += tmpNode[j]/NP;
      }
    }
  }

  int patch(int i, int j)
  {
    auto & loc = m_top.locations();
    auto & nei = m_top.neighbors();
    return nei[loc[i]+j];
  }

  int patch_size(int i)
  {
    auto & loc = m_top.locations();
    return loc[i+1] -loc[i];
  }
  
private:
  Toplogy m_top;
  std::shared_ptr<PMesh> m_mesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC

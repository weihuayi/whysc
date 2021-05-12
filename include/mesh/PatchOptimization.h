#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <numeric>

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg;

template<typename PMesh, typename MeshQuality>
class PatchOptimization
{
public:
  typedef typename PMesh::Node Node;
  typedef typename PMesh::Vector Vector;
  typedef typename PMesh::Toplogy Toplogy;

public:
  PatchOptimization(std::shared_ptr<PMesh> mesh)
  {
    m_mesh = mesh;
    m_quality = std::make_shared<MeshQuality>(mesh);
    mesh->node_to_cell(m_top);
  }

  void optimization(int i, Node & node, std::string method = "bar")
  {
    if(method=="bar")
      bar_optimization(i, node);
    else if(method=="tet")
      tet_optimization(i, node);
  }

  void bar_optimization(int i, Node & node)
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

  void tet_optimization(int i, Node & node)
  {
    int NP = patch_size(i);
    std::vector<int> cidx(NP);
    Vector v = {0};
    for(int k = 0; k < NP; k++)
    {
      cidx[k] = patch(i, k);
      Vector tmpVector;
      double w;
      m_quality->nabla_quality(patch(i, k), local_index(i, k), tmpVector, w);
      v = v + w*tmpVector;
    }

    Node tmpNode(m_mesh->nodes(i));

    double a = 0;
    double b = 10;//TODO 把这个改为边长
    double Qa = m_quality->patch_quality(cidx);
    m_mesh->nodes(i) = tmpNode + b*v;
    double Qb = m_quality->patch_quality(cidx);

    while(abs(b-a)>0.001)//TODO 0.618法
    {
      int c = (a+b)/2;
      m_mesh->nodes(i) = tmpNode + c*v;
      double Qc = m_quality->patch_quality(cidx);

      L*v;
      m_mesh
    }
  }

  int local_index(int i, int j)
  {
    auto & loc = m_top.locations();
    auto & locidx = m_top.local_indices();
    return locidx[loc[i]+j];
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
  std::shared_ptr<MeshQuality> m_quality;
};

} // end of namespace Mesh

} // end of namespace WHYSC

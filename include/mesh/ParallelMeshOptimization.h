#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include "PatchOptimization.h"

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg;

template<typename PMesh, typename MeshQuality>
class ParallelMeshOptimization
{
public:
  typedef typename PMesh::Node Node;
  typedef typename PMesh::NodeArray NodeArray;
  typedef typename PMesh::Toplogy Toplogy;
  typedef GhostFillingAlg<PMesh> SetGhostAlg;
  typedef PatchOptimization<PMesh, MeshQuality> PatchOpt;

public:
  ParallelMeshOptimization(std::shared_ptr<PMesh> mesh, std::vector<int> & color, int cmax, MPI_Comm comm)
  {
    m_cmax = cmax;
    m_comm = comm;
    m_mesh = mesh;
    m_color = color;
    m_node = mesh->nodes();
    m_set_ghost_alg = std::make_shared<SetGhostAlg>(mesh, comm);
    m_patch = std::make_shared<PatchOpt>(mesh);
  }

  void mesh_optimization(std::string method="bar")
  {
    auto GD = m_mesh->geo_dimension();
    auto NN = m_mesh->number_of_nodes();
    auto NC = m_mesh->number_of_cells();
    auto gdof = m_mesh->nodedata().gdof;
    auto & isghostnode = m_set_ghost_alg->get_ghost_node();

    //int cMax;
    //auto cmax = *max_element(m_color.begin(), m_color.end());
    //MPI_Allreduce(&cmax, &cMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for(int i = 1; i <= m_cmax; i++)
    {
      for(int j = 0; j < NN; j++)
      {
        if((m_color[j] == i) & (gdof[j] > 0) & (!isghostnode[j]))
        {
          m_patch->optimization(j, m_node[j], method);
        }
      }
      m_set_ghost_alg->fill(m_node, GD);
      m_mesh->nodes() = m_node;
    }
  }
 
private:
  int m_cmax;
  MPI_Comm m_comm;
  NodeArray m_node;
  std::shared_ptr<PMesh> m_mesh;
  std::shared_ptr<PMesh> m_auxmesh;
  std::vector<int> m_color;
  std::shared_ptr<GhostFillingAlg<PMesh> > m_set_ghost_alg;
  std::shared_ptr<PatchOpt> m_patch;
};

} // end of namespace Mesh

} // end of namespace WHYSC

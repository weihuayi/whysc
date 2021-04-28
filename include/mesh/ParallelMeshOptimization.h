#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg;

template<typename PMesh>
class ParallelMeshOptimization
{
public:
  typedef typename PMesh::Node Node;
  typedef typename PMesh::Toplogy Toplogy;
  typedef GhostFillingAlg<PMesh> SetGhostAlg;

public:
  ParallelMeshOptimization(std::shared_ptr<PMesh> mesh, std::vector<int> & color, MPI_Comm comm)
  {
    m_comm = comm;
    m_mesh = mesh;
    m_color = color;
    m_node = mesh->nodes();
    m_set_ghost_alg = std::make_shared<SetGhostAlg>(mesh, comm);
  }

  void mesh_optimization()
  {
    Toplogy node2cell;
    m_mesh->node_to_cell(node2cell);
    auto & loc = node2cell.locations();
    auto & nei = node2cell.neighbors();

    auto GD = m_mesh->geo_dimension();
    auto NN = m_mesh->number_of_nodes();
    auto NC = m_mesh->number_of_cells();
    auto fixed = m_mesh->data()["fixednode"];

    int cMax;
    auto cmax = *max_element(m_color.begin(), m_color.end());
    MPI_Allreduce(&cmax, &cMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for(int i = 1; i <= cMax; i++)
    {
      for(int j = 0; j < NN; j++)
      {
        if(m_color[j] == i & fixed[j] == 0)
        {
          int N = loc[j+1] - loc[j];
          std::vector<int> patch(N);
          for(int k = 0; k < N; k++)
          {
            patch[k] = nei[loc[j]+k];
          }
          patch_optimization(j, patch);
        }
      }
      m_set_ghost_alg->fill(m_node, GD);
      m_mesh->nodes() = m_node;
    }
  }

  void patch_optimization(int i, std::vector<int> & patch)
  {
    int GD = m_mesh->geo_dimension();
    int NP = patch.size();
    Node node{0, 0, 0}; 
    for(int k = 0; k < NP; k++)
    {
      Node tmpNode;
      m_mesh->cell_barycenter(path[k], tmpNode);
      for(int j = 0; j < GD; j++)
      {
        node[j] += tmpNode[j]/NP;
      }
    }
    m_node[i] = node;
  }
  
private:
  MPI_Comm m_comm;
  std::shared_ptr<PMesh> m_mesh;
  std::vector<int> m_color;
  std::vector<Node> m_node;
  std::shared_ptr<GhostFillingAlg<PMesh> > m_set_ghost_alg;
};

} // end of namespace Mesh

} // end of namespace WHYSC

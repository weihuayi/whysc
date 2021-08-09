#ifndef ParallelMeshOptAlg_h
#define ParallelMeshOptAlg_h

#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>

#include "NodePatchOptAlg.h"
#include "ParallelMeshColoringAlg.h"

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg;

template<typename PMesh, typename ObjectFunction, typename Model>
class ParallelMeshOptAlg
{
public:
  typedef typename PMesh::Node Node;
  typedef typename PMesh::NodeArray NodeArray;
  typedef typename PMesh::Toplogy Toplogy;
  typedef GhostFillingAlg<PMesh> SetGhostAlg;
  typedef NodePatchOptAlg<PMesh, ObjectFunction, Model> PatchOpt;
  typedef ParallelMeshColoringAlg<PMesh> PMeshColoring;

public:
  ParallelMeshOptAlg(std::shared_ptr<PMesh> mesh, std::shared_ptr<Model> model, 
      MPI_Comm comm): m_comm(comm), m_mesh(mesh)
  {
    m_set_ghost_alg = std::make_shared<SetGhostAlg>(mesh, comm);
    m_coloring_alg = std::make_shared<PMeshColoring>(mesh, m_set_ghost_alg, comm);
    m_node_patch_opt_alg = std::make_shared<PatchOpt>(mesh, model);
    m_coloring_alg->coloring(true);
  }

  void optimization()
  {
    auto GD = m_mesh->geo_dimension();
    auto NN = m_mesh->number_of_nodes();
    auto & gdof = m_mesh->get_node_int_data()["gdof"];
    auto & color2node = m_coloring_alg->get_color_to_node();

    for(auto & idxs : color2node)//循环所有的颜色
    {
      for(auto & i : idxs)//循环所有的点
      {
        if(gdof[i] > 0)//只优化当前颜色自由度大于0, 且是本网格的点
        {

          m_node_patch_opt_alg->optimization(i);//优化节点
        }
      }
      m_set_ghost_alg->fill(m_mesh->nodes(), GD); //通信重叠区节点位置
    }
  }

private:
  MPI_Comm m_comm;
  std::shared_ptr<PMesh> m_mesh;
  std::shared_ptr<SetGhostAlg> m_set_ghost_alg;
  std::shared_ptr<PMeshColoring> m_coloring_alg;
  std::shared_ptr<PatchOpt> m_node_patch_opt_alg;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMeshOptAlg_h

#include <memory>
#include <mpi.h>
#include <vector>

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class Communication 
{
public:
  typedef typename PMesh::I I;
  typedef typename PMesh::Toplogy Toplogy;

public:
  Communication(std::shared_ptr<PMesh> mesh, MPI_Comm comm)
  {
    int NN = mesh->number_of_nodes();
    m_comm = comm;
    m_mesh = mesh;

    auto & npid = mesh->node_process_id();
    auto & isGhostNode = get_ghost_node();
    isGhostNode.resize(NN);
    for(int i = 0; i < NN; i++)
    {
      if(npid[i]==mesh->id())
        isGhostNode[i] = false;
      else
        isGhostNode[i] = true;
    }
  }

  template<typename F>
  void communicate(std::vector<F> & data)
  {
    auto mesh = get_mesh();
    auto & isGhostNode = get_ghost_node();
    auto & pds = mesh->parallel_data_structure();
    for(auto map : pds)
    {
      auto & target = map.first; 
      auto & meshOverlap = map.second; 
      auto & overlap = meshOverlap.entity_overlap(0);

      auto & locid = overlap.loc_index();
      auto & adjid = overlap.adj_index();

      int N = adjid.size();
      int adjData[2*N];

      for(int j = 0; j < N; j++)
      {
        adjData[j*2] = -1;
        if(!isGhostNode[locid[j]])//只发送自己的数据
        {
          adjData[j*2] = adjid[j];
          adjData[j*2+1] = data[locid[j]];
        }
      }
      MPI_Send(adjData, N*2, MPI_INT, target, 1, m_comm);
    }//发送数据完成

    for(auto map : pds)
    {
      auto & target = map.first; 
      auto & meshOverlap = map.second; 
      auto & overlap = meshOverlap.entity_overlap(0);
      auto & adjid = overlap.adj_index();

      int N = adjid.size();
      int locData[2*N];

      MPI_Recv(locData, N*2, MPI_INT, target, 1, m_comm, MPI_STATUS_IGNORE);
      for(int k = 0; k < N; k++)
      {
        if(locData[2*k] >= 0)
        {
          if(isGhostNode[locData[k*2]])//只接收别人的数据
          {
            data[locData[2*k]] = locData[k*2+1];//填充影像节点数据
          }
        }
      }
    }//接收数据完成
  }

  std::vector<bool> & get_ghost_node()
  {
    return m_isGhostNode;
  }

  std::shared_ptr<PMesh> get_mesh()
  {
    return m_mesh;
  }

private:
  MPI_Comm m_comm;
  std::shared_ptr<PMesh> m_mesh;
  std::vector<bool> m_isGhostNode;
};
} // end of namespace Mesh

} // end of namespace WHYSC

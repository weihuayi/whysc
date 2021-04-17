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

    auto & pds = mesh->parallel_data_structure();
    auto & isImageNode = get_image_node();
    std::vector<bool> isBdNode;
    mesh->is_Boundary_node(isBdNode);

    isImageNode.resize(NN);
    std::vector<bool> isOverlapNode(NN, 0);

    for(auto & map : pds)
    {
      auto & meshOverlap = map.second; 
      auto & overlap = meshOverlap.entity_overlap(0);

      auto & locid = overlap.loc_index();
      for(auto id : locid)
      {
        isOverlapNode[id] = true;
      }
    }

    Toplogy node2node;
    mesh->node_to_node(node2node);
    auto & loc = node2node.locations();
    auto & nei = node2node.neighbors();
    for(int i = 0; i < NN; i++)
    {
      isImageNode[i] = false;
      if(isOverlapNode[i] & isBdNode[i])
      {
        isImageNode[i] = true;
        for(int j = loc[i]; j < loc[i+1]; j++)
        {
          isImageNode[i] = isImageNode[i] & isOverlapNode[nei[j]];
        }
      }
    }
  }

  template<typename F>
  void communicate(std::vector<F> & data)
  {
    auto mesh = get_mesh();
    auto & isImageNode = get_image_node();
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
        if(!isImageNode[locid[j]])//只发送自己的数据
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
          if(isImageNode[locData[k*2]])//只接收别人的数据
          {
            data[locData[2*k]] = locData[k*2+1];//填充影像节点数据
          }
        }
      }
    }//接收数据完成
  }

  std::vector<bool> & get_image_node()
  {
    return m_isImageNode;
  }

  std::shared_ptr<PMesh> get_mesh()
  {
    return m_mesh;
  }

private:
  MPI_Comm m_comm;
  std::shared_ptr<PMesh> m_mesh;
  std::vector<bool> m_isImageNode;
};
} // end of namespace Mesh

} // end of namespace WHYSC

#ifndef GhostFillingAlg_h
#define GhostFillingAlg_h

#include <memory>
#include <iostream>
#include <mpi.h>
#include <map>
#include <vector>
#include <array>

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class GhostFillingAlg 
{
public:
  typedef typename PMesh::I I;
  typedef typename PMesh::Toplogy Toplogy;

public:
  GhostFillingAlg(std::shared_ptr<PMesh> mesh, MPI_Comm comm)
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

    auto & pds = mesh->parallel_data_structure();
    for(auto & map : pds)
    {
      auto & target = map.first; 
      auto & meshOverlap = map.second; 
      auto & overlap = meshOverlap.entity_overlap(0);

      auto & locid = overlap.loc_index();
      auto & adjid = overlap.adj_index();

      m_dataid[target][0] = &locid;
      m_dataid[target][1] = &adjid;
    }
  }

  template<typename F>
  void fill(std::vector<F> & data, int n)
  {
    int NP;
    MPI_Comm_size(MPI_COMM_WORLD, &NP);

    auto mesh = get_mesh();
    auto & isGhostNode = get_ghost_node();

    for(auto & map : m_dataid)
    {
      auto & target = map.first; 
      auto & locid = map.second[0];
      auto & adjid = map.second[1];

      int N = adjid->size();
      double adjData[(n+1)*N];
      double locData[(n+1)*N];
      for(int j = 0; j < N; j++)
      {
        adjData[j*(n+1)] = -1;
        if(!isGhostNode[locid->at(j)])//只发送自己的数据
        {
          adjData[j*(n+1)] = (double)adjid->at(j);
          for(int k = 1; k < n+1; k++)
            adjData[j*(n+1)+k] = (double)data[locid->at(j)][k-1];
        }
      }
      //发送并接受 target 的数据
      MPI_Sendrecv(adjData, N*(n+1), MPI_DOUBLE, target, 1, 
          locData, N*(n+1), MPI_DOUBLE, target, 1, m_comm, MPI_STATUS_IGNORE);

      for(int k = 0; k < N; k++)
      {
        if(locData[(n+1)*k] >= 0)
        {
          if(isGhostNode[(int)locData[k*(n+1)]])//只接收别人的数据
          {
            for(int j = 1; j < n+1; j++)
              data[(int)locData[(n+1)*k]][j-1] = locData[k*(n+1)+j];//填充影像节点数据
          }
        }
      }
    }
  }

  template<typename F>
  void fill(std::vector<F> & data)
  {
    int NP;
    MPI_Comm_size(MPI_COMM_WORLD, &NP);

    auto mesh = get_mesh();
    auto & isGhostNode = get_ghost_node();

    for(auto & map : m_dataid)
    {
      auto & target = map.first; 
      auto & locid = map.second[0];
      auto & adjid = map.second[1];

      int N = adjid->size();
      double adjData[2*N];
      double locData[2*N];
      for(int j = 0; j < N; j++)
      {
        adjData[j*2] = -1;
        if(!isGhostNode[locid->at(j)])//只发送自己的数据
        {
          adjData[j*2] = (double)adjid->at(j);
          adjData[j*2+1] = (double)data[locid->at(j)];
        }
      }
      MPI_Sendrecv(adjData, N*2, MPI_DOUBLE, target, 1, 
          locData, N*2, MPI_DOUBLE, target, 1, m_comm, MPI_STATUS_IGNORE);

      for(int k = 0; k < N; k++)
      {
        if(locData[2*k] >= 0)
        {
          if(isGhostNode[(int)locData[k*2]])//只接收别人的数据
          {
            data[(int)locData[2*k]] = locData[k*2+1];//填充影像节点数据
          }
        }
      }
    }//发送数据完成

    /*
    for(auto & map : m_dataid)
    {
      auto & target = map.first; 
      auto & adjid = map.second[1];

      int N = adjid->size();
      double locData[2*N];

      std::cout<< mesh->id() << " recving " << sizeof(locData) << " to " <<target<<std::endl; 
      MPI_Recv(locData, N*2, MPI_DOUBLE, target, 1, m_comm, MPI_STATUS_IGNORE);
      std::cout<< mesh->id() << " recved " << sizeof(locData) << " to " <<target<<std::endl; 
      for(int k = 0; k < N; k++)
      {
        if(locData[2*k] >= 0)
        {
          if(isGhostNode[(int)locData[k*2]])//只接收别人的数据
          {
            data[(int)locData[2*k]] = locData[k*2+1];//填充影像节点数据
          }
        }
      }
    }//接收数据完成
    */
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
  std::map<int, std::array<std::vector<int>*, 2> > m_dataid;
  std::vector<bool> m_isGhostNode;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of GhostFillingAlg_h

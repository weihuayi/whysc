#ifndef ParallelMesher_h
#define ParallelMesher_h

#include <string>
#include <memory>
#include <mpi.h>

#include "VTKMeshReader.h"

namespace WHYSC {
namespace Mesh {

template<typename PMesh>
class ParallelMesher 
{
public:
  typedef VTKMeshReader<PMesh> Reader;
  typedef typename PMesh::Toplogy Toplogy;
  typedef typename PMesh::I I;

public:
  ParallelMesher(std::string fnamebase, std::string fnameextent, MPI_Comm comm = MPI_COMM_WORLD)
  {
    m_comm = comm;
    MPI_Comm_rank(m_comm, &m_rank);
    m_pmesh = std::make_shared<PMesh>(m_rank);


    std::stringstream ss;
    ss << fnamebase << "_" << m_rank << fnameextent;
    std::string fname = ss.str();

    Reader reader(m_pmesh);
    reader.read(ss.str());

    auto & gid = m_pmesh->node_global_id();
    reader.get_node_data("gid", gid);

    std::vector<int> npid(gid.size());
    reader.get_node_data("nid", npid);

    build_mesh(npid);
  }

  //virtual ~ParallelMesher();

  void build_mesh(std::vector<int> & npid)
  {

    Toplogy node2node;
    m_pmesh->node_to_node(node2node);
    auto & loc = node2node.locations();
    auto & nei = node2node.neighbors();

    auto id = m_pmesh->id();
    auto NN = m_pmesh->number_of_nodes();
    auto & pds = m_pmesh->parallel_data_structure();
    auto & gid = m_pmesh->node_global_id();
    auto & lnn = m_pmesh->number_of_local_nodes();
    lnn = 0;

    std::map<I, std::map<I, I> > ng2l; //每个重叠区的节点全局编号到局部编号的映射
    for(I i=0; i < NN; i++)
    {
      //不是本进程的点的周围的点就是 overlap 的点
      if(npid[i] != id) // i 是其他进程的点, 那么他相邻的点也是 npid[i] 网格的点
      {
        ng2l[npid[i]][gid[i]] = i;
        for(int k = loc[i]; k < loc[i+1]; k++)
        {
          ng2l[npid[i]][gid[nei[k]]] = nei[k];
        }
      }
      else
      {
        lnn++;
      }
    }//完成 ng2l 的构建

    //发送 ng2l 中全局与局部的编号信息
    for(auto map : ng2l)
    {
      auto target  = map.first;
      auto & idxmap = map.second;

      int N = idxmap.size();
      int data[N*2];
      int j = 0;
      for(auto pair : idxmap)
      {
        data[2*j] = pair.first;
        data[2*j+1] = pair.second; //传输的数据是先全局后局部
        j++;
      }

      MPI_Send(data, N*2, MPI_INT, target, 1, m_comm);
    }

    //接收重叠区的编号信息
    for(auto map : ng2l)
    {
      auto target  = map.first;
      auto & idxmap = map.second;

      int N = idxmap.size();
      int data[N*2];

      MPI_Recv(data, N*2, MPI_INT, target, 1, m_comm, MPI_STATUS_IGNORE);

      //把 data 中的数据和 idxmap 结合
      pds[target].init(3);
      auto & overlap0 = pds[target].entity_overlap(0);
      auto & locid = overlap0.loc_index();
      auto & adjid = overlap0.adj_index();

      locid.resize(N);
      adjid.resize(N);

      for(int j = 0; j < N; j++)
      {
        locid[j] = idxmap[data[2*j]];
        adjid[j] = data[2*j+1];
      }
    }
  }

  std::shared_ptr<PMesh> get_mesh()
  {
    return m_pmesh;
  }

private:

  MPI_Comm m_comm;
  int m_rank;
  std::shared_ptr<PMesh> m_pmesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of ParallelMesher_h

#ifndef ParallelMesher_h
#define ParallelMesher_h

#include <string>
#include <memory>
#include <array>
#include <vector>
#include <mpi.h>
#include <map>
#include <set>

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

    auto & npid = m_pmesh->node_process_id();
    reader.get_node_data("nid", npid);

    auto & nodeIntData = m_pmesh->get_node_int_data();
    reader.get_node_data("gdof", nodeIntData["gdof"]);
    reader.get_node_data("gtag", nodeIntData["gtag"]);
    nodeIntData["nid"] = npid;
    build_mesh();
  }

  //virtual ~ParallelMesher();

  void build_mesh()
  {
    auto mesh = get_mesh();
    auto & cell = mesh->cells();

    auto id = mesh->id();
    auto NN = mesh->number_of_nodes();
    auto & pds = mesh->parallel_data_structure();
    auto & gid = mesh->node_global_id();
    auto & npid = mesh->node_process_id();
    std::map<I, std::map<I, I> > ng2l; //每个重叠区的节点全局编号到局部编号的映射
    for(auto & c : cell)
    {
      //不是本进程的点的周围的点就是 overlap 的点
      std::set<int> idx;
      for(auto v: c)
      {
        if(npid[v]!=id)
          idx.insert(npid[v]); // cell 在其节点所在网格
      }

      for(auto i: idx)
      {
        for(auto v : c)
        {
          ng2l[i][gid[v]] = v;
        }
      }
    }//完成 ng2l 的构建

    communicate_index(ng2l, 0);
    //build_overlap(mesh->edges(), 1);
    //build_overlap(mesh->cells(), 2);
  }

  void communicate_index(std::map<I, std::map<I, I> > & g2l, int d)
  {
    auto mesh = get_mesh();
    auto & pds = mesh->parallel_data_structure();

    for(auto map : g2l)
    {
      auto target  = map.first;
      auto & idxmap = map.second;

      int N = idxmap.size();
      int mdata[N*2];
      int odata[N*2];
      int j = 0;
      for(auto pair : idxmap)
      {
        mdata[2*j] = pair.first;
        mdata[2*j+1] = pair.second; //传输的数据是先全局后局部
        j++;
      }

      MPI_Sendrecv(mdata, N*2, MPI_INT, target, 1, 
          odata, N*2, MPI_INT, target, 1, m_comm, MPI_STATUS_IGNORE);

      pds[target].init(3);
      auto & overlap0 = pds[target].entity_overlap(d);
      auto & locid = overlap0.loc_index();
      auto & adjid = overlap0.adj_index();

      locid.resize(N);
      adjid.resize(N);

      for(int j = 0; j < N; j++)
      {
        locid[j] = idxmap[odata[2*j]];
        adjid[j] = odata[2*j+1];
      }
    }
  }

  template<typename Vec>
  void build_overlap(Vec & com, int d)
  {
    auto mesh = get_mesh();
    int NN = mesh->number_of_nodes();

    auto & pds = mesh->parallel_data_structure();
    std::map<I, std::map<I, I> > g2l; //每个重叠区的 d 维单形全局编号到局部编号的映射

    for(auto & map : pds)
    {
      auto & target = map.first;
      auto & meshOverlap = map.second; 
      auto & overlap0 = meshOverlap.entity_overlap(0);
      auto & overlapd = meshOverlap.entity_overlap(d);

      std::vector<bool> isOverlapNode(NN, 0);
      auto & locid = overlap0.loc_index();
      for(auto id : locid)
      {
        isOverlapNode[id] = true;
      }

      for(int i = 0; i < com.size(); i++)
      {
        bool flag = true;

        int E[d+1];
        for(int j = 0; j < d+1; j++)
        {
          flag = flag & isOverlapNode[com[i][j]]; //顶点都在 overlap 则自己也在 overlap
          E[j] = com[i][j];
        }
        std::sort(E, E+d+1);

        if(flag)// i 单形在 overlap 中
        {
          int gid = 0;
          for(int k = 1; k < d+2; k++)
          {
            int num0 = 1;
            int num1 = 1;
            for(int j = 0; j < k; j++)
            {
              num0 *= E[k-1]+j;
              num1 *= j+1;
            }
            gid += num0/num1;
          }
          g2l[target][gid] = i;
        }
      }
    }// g2l 构建完成
    communicate_index(g2l, d);
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

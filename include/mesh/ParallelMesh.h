#ifndef ParallelMesh_h
#define ParallelMesh_h

#include "thirdparty/json.hpp"

namespace WHYSC {
namespace Mesh {

template<typename GK, typename Mesh>
class ParallelMesh: public Mesh
{
public:
  typedef typename GK::Int I;
  typedef typename GK::Float F;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Edge Edge;
  typedef typename Mesh::Face Face;
  typedef typename Mesh::Cell Cell;

  typedef typename Mesh::Face2cell Face2cell;
  typedef typename Mesh::Cell2face Cell2cell;
  typedef typename Mesh::Cell2cell Cell2cell;
  typedef typename Mesh::Cell2edge Cell2edge;

  typedef typename Mesh::Toplogy Toplogy;

  typedef typename Mesh::NodeIterator NodeIterator;
  typedef typename Mesh::EdgeIterator EdgeIterator;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::CellIterator CellIterator;

  typedef std::map<I, std::set<I> > PDS; // 并行网格数据结构
  typedef json MeshInfo;

public:
  ParallelMesh(int id, std::string & ptype="node")
  {
    m_info["id"] = id;
    m_info["ptype"] = ptype;

    m_info["cpid"] = std::vector<I>(); // 单元所在进程编号
    m_info["npid"] = std::vector<I>(); // 节点所在进程编号

    m_info["cgid"] = std::vector<I>(); // 单元的全局编号
    m_info["ngid"] = std::vector<I>(); // 节点的全局编号

    m_info["cg2l"] = std::map<I, I>(); // 单元全局到局部编号的映射
    m_info["ng2l"] = std::map<I, I>(); // 节点全局到局部编号的映射

    m_info["pds"]["node"] = PDS();
    m_info["pds"]["cell"] = PDS();
  }

  void construct_parallel_data_structure()
  {

    auto & ptype = m_info["ptype"];
    auto id = m_info["id"];
    if(ptype == "node")
    {
      Toplogy node2node;
      node_to_node(node2node);
      auto & loc = node2node.locations();
      auto & nei = node2node.neighbors();
      auto & pds = parallel_data_structure();// 本网格需要发送信息的节点 
      auto NN = number_of_nodes();
      auto & npid = node_process_id();
      auto & gid = node_global_id();
      auto & ng2l = node_global_to_local_id();
      for(I i=0; i < NN; i++)
      {
        ng2l[gid[i]] = i;
        if(npid[i] != id) // i 是其他进程的点, 那么他相邻的点也是 nid[j] 网格的点
        {
          for(int k = loc[i]; k < loc[i+1]; k++)// 循环 i 相邻的点, 将本进程的点放入pds0[nid[j]] 中
          {
            if(npid[k] == id) // k 是本进程的点
            {
              pds[npid[i]].insert(k);
            }
          }
        }
      }
    }
    else if(ptype == "cell")
    {
    }
  }

  std::map<I, I> & cell_global_to_local_id()
  {
    return m_info["cg2l"];
  }

  std::map<I, I> & node_global_to_local_id()
  {
    return m_info["ng2l"];
  }

  std::vector<I> & cell_global_id()
  {
    return m_info["cgid"]
  }
  
  std::vector<I> & node_global_id()
  {
    return m_info["ngid"];
  }

  std::vector<I> & node_process_id()
  {
    return m_info["npid"];
  }

  std::vector<I> & cell_process_id()
  {
    return m_info["cpid"];
  }

  PDS & parallel_data_structure()
  {
    auto & ptype = m_info["ptype"];
    return m_info["pds"][ptype];
  }
private:
  MeshInfo m_info;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

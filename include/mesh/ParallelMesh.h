#ifndef ParallelMesh_h
#define ParallelMesh_h

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
  typedef typename Mesh::Cell2face Cell2face;
  typedef typename Mesh::Cell2edge Cell2edge;

  typedef typename Mesh::Toplogy Toplogy;

  typedef typename Mesh::NodeIterator NodeIterator;
  typedef typename Mesh::EdgeIterator EdgeIterator;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::CellIterator CellIterator;

  typedef std::map<I, std::set<I> > PDS; // 并行网格数据结构

public:
  ParallelMesh(int id, int gw=0)
  {
    m_id = id;
    m_gw = gw; 
  }

  void construct_parallel_data_structure()
  {

    auto id = m_id;

    Toplogy node2node;
    this->node_to_node(node2node);
    auto & loc = node2node.locations();
    auto & nei = node2node.neighbors();
    auto & pds = parallel_data_structure();// 本网格需要发送信息的节点 
    auto NN = this->number_of_nodes();
    auto & npid = node_process_id();
    auto & gid = node_global_id();
    auto & ng2l = node_global_to_local_id();
    //auto & num = number_of_nodes_in_process();
    for(I i=0; i < NN; i++)
    {
      //num[npid[i]]++;
      ng2l[gid[i]] = i;
      if(npid[i] != id) // i 是其他进程的点, 那么他相邻的点也是 nid[j] 网格的点
      {
        for(int k = loc[i]; k < loc[i+1]; k++)//循环i相邻的点, 将本进程的点放入pds0[nid[j]] 中
        {
          if(npid[k] == id) // k 是本进程的点
          {
            pds[npid[i]].insert(k);
          }
        }
      }
    }
  }

  int id()
  {
    return m_id;
  }

  std::vector<I> & number_of_nodes_in_process()
  {
    return m_NN;
  }

  std::vector<I> & number_of_cells_in_process()
  {
    return m_NC;
  }

  std::map<I, I> & cell_global_to_local_id()
  {
    return m_cg2l;
  }

  std::map<I, I> & node_global_to_local_id()
  {
    return m_ng2l;
  }

  std::vector<I> & cell_global_id()
  {
    return m_cgid;
  }
  
  std::vector<I> & node_global_id()
  {
    return m_ngid;
  }

  std::vector<I> & node_process_id()
  {
    return m_npid;
  }

  std::vector<I> & cell_process_id()
  {
    return m_cpid;
  }

  PDS & parallel_data_structure()
  {
    return m_pds;
  }

private:
  I m_LNN;
  I m_LNC;
  I m_id; // 网格块编号
  int m_gw; // 影像区宽度
  std::vector<I> m_NN;
  std::vector<I> m_NC;
  std::vector<I> m_cpid; //单元所在进程编号  
  std::vector<I> m_npid; //节点所在进程编号
  std::vector<I> m_cgid; //单元的全局编号
  std::vector<I> m_ngid; //节点的全局编号
  std::map<I, I> m_cg2l; //单元全局到局部编号的映射
  std::map<I, I> m_ng2l; //节点全局到局部编号的映射

  PDS m_pds;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

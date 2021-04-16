#ifndef ParallelMesh_h
#define ParallelMesh_h

#include <vector>
#include <map>

namespace WHYSC {
namespace Mesh {

template<typename I>
class EntityOverlap
{
public:
  typedef std::vector<I> Container;

public:
  EntityOverlap()
  {
    m_init_flag = false;
  }

  bool empty()
  {
    return m_init_flag;
  }

  // 实体集在本进程的编号数组
  Container & loc_index()
  {
    return m_loc;
  }

  // 实体集在邻接进程的编号数组
  Container & adj_index()
  {
    return m_adj;
  }

private:
  bool m_init_flag; // 是否已经初始化
  Container m_loc; // 实体在当前网格块的编号
  Container m_adj; // 实体在邻居网格块的编号
};

template<typename I>
class MeshOverlap
{
public:

  MeshOverlap()
  {
    m_GD = -1;
  }

  MeshOverlap(int GD)
  {
    init(GD);
  }

  void init(int GD)
  {
    m_GD = GD;
    m_entity_overlap.resize(GD);
  }

  bool empty()
  {
    return m_GD == -1;
  }

  I get_geo_dimension() const
  {// 网格重叠的最高维度
    return m_GD;
  }

  void set_geo_dimension(I GD)
  {// 网格重叠的最高维度
    m_GD = GD;
  }

  // 返回第 i 维的重叠实体编号
  EntityOverlap<I> & entity_overlap(I i)
  {
    // 0 <= i <= m_GD
    return m_entity_overlap[i];
  }

private:
  I m_GD; // 相邻网格块重叠实体的最高维度 
  std::vector<EntityOverlap<I> > m_entity_overlap; // 长度是 m_GD
};

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

  // 并行网格块之间的拓扑关系, 记录了当前网格块和相邻网格块之间重叠网格实体
  // 在本网格块和相邻网格块之间的对应关系
  typedef std::map<I, MeshOverlap<I> > PDS; 

public:

  ParallelMesh(int id)
  {
    m_id = id;
  }

  void construct_parallel_data_structure(std::vector<int> & npid, MPI_Comm & comm = MPI_COMM_WORLD)
  {
    auto id = m_id;

    Toplogy node2node;
    Mesh::node_to_node(node2node);
    auto & loc = node2node.locations();
    auto & nei = node2node.neighbors();

    auto NN = Mesh::number_of_nodes();
    auto & pds = parallel_data_structure();
    auto & gid = node_global_id();
    auto & lnn = number_of_local_nodes();
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

      MPI_Send(data, N*2, MPI_INT, target, 1, comm);
      //std::cout<< m_id << "已经发送给" << target << "信息" <<std::endl;
    }

    //接收重叠区的编号信息
    for(auto map : ng2l)
    {
      auto target  = map.first;
      auto & idxmap = map.second;

      int N = idxmap.size();
      std::cout<< "N = " << N <<std::endl;
      int data[N*2];

      MPI_Recv(data, N*2, MPI_INT, target, 1, comm, MPI_STATUS_IGNORE);
      //std::cout<< m_id << "接收到" << target << "的信息" <<std::endl;

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

    if(m_id==0)
    {
      auto & overlap0 = pds[1].entity_overlap(0);
      auto & locid = overlap0.loc_index();
      auto & adjid = overlap0.adj_index();
      std::cout<< "locid" <<std::endl;
      for(auto a : locid)
        std::cout<< a <<std::endl;

      std::cout<< "adjid" <<std::endl;
      for(auto a : adjid)
        std::cout<< a <<std::endl;
    }


  }

  int id()
  {
    return m_id;
  }

  int & number_of_local_nodes()
  {
    return m_lnn;
  }

  std::vector<I> & cell_global_id()
  {
    return m_cgid;
  }
  
  std::vector<I> & node_global_id()
  {
    return m_ngid;
  }

  PDS & parallel_data_structure()
  {
    return m_pds;
  }

private:
  I m_id; // 网格块编号(进程的编号）
  I m_lnn;
  std::vector<I> m_cgid; //单元的全局编号
  std::vector<I> m_ngid; //节点的全局编号
  PDS m_pds;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

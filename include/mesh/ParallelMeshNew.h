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

  ParallelMesh(int id): m_id(id) {}

  ParallelMesh(): m_id(0) {}

  int id()
  {
    return m_id;
  }

  std::vector<I> & node_process_id()
  {
    return m_npid;
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
  std::vector<I> m_cgid; //单元的全局编号
  std::vector<I> m_ngid; //节点的全局编号
  std::vector<I> m_npid; //节点的进程编号
  PDS m_pds;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

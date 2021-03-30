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
  typedef json MeshInfo;
public:
  ParallelMesh(int id, std::string & name, std::string & ptype="node")
  {
    m_info["name"] = name;
    m_info["ptype"] = ptype;
    m_info["id"] = id;
  }

  std::vector<I> & cell_global_id()
  {
  }
  
  std::vector<I> & node_global_id()
  {
    return m_node_gid;
  }

  ParallelDataStructure & parallel_data_structure()
  {
    return m_pds;
  }
private:
  MeshInfo m_info;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

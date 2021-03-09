#ifndef ParallelMesh_h
#define ParallelMesh_h

#include <mpi.h>

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

public:
  ParallelMesh(MPI_Comm comm=MPI_COMM_WORLD)
  {
    m_comm = comm;
    MPI_Comm_size(m_comm, &m_wsize);
    MPI_Comm_rank(m_comm, &m_wrank);
  }
  
private:
  std::vector<I> m_data; 
  int m_wsize; // world size
  int m_wrank; // world rank
  MPI_Comm m_comm;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ParallelMesh_h

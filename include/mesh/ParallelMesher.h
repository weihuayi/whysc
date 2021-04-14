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
public:
  ParallelMesher(MPI_Comm comm=MPI_COMM_WORLD, std::string fnamebase, std::string fnameextent)
  {
    m_comm = comm;
    MPI_Comm_rank(m_comm, &m_rank);
    m_pmesh = std::make_shared<PMesh>(m_rank);

    std::stringstream ss;
    ss << fnamebase << "_" << rank << fnameextent;
    std::string fname = ss.str();
  Reader reader(&mesh);
  reader.read(ss.str());

  }

  virtual ~ParallelMesher();

  std::shared_ptr<Pmesh> build_mesh()
  {
  }
  
private:
  Reader m_reader;
  MPI_Comm m_comm;
  int m_rank;
  std::shared_ptr<PMesh> m_pmesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of ParallelMesher_h

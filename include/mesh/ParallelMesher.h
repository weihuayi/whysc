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

    m_pmesh->construct_parallel_data_structure(npid, m_comm);
  }

  //virtual ~ParallelMesher();

  std::shared_ptr<PMesh> build_mesh()
  {
    return m_pmesh;
  }
  
private:
  //Reader m_reader;
  MPI_Comm m_comm;
  int m_rank;
  std::shared_ptr<PMesh> m_pmesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of ParallelMesher_h

#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/VTKMeshReader.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PMesh;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;


void mesh_coloring(TriMesh & mesh, std::vector<int> & color)
{
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto NN = mesh.number_of_nodes();
  auto NE = mesh.number_of_edges();

  color.resize(NN);

  std::vector<int> randVal(NN);
  std::vector<bool> isMin(NN, true);

  std::list<int> edges; //没有被删除的边
  std::list<int> nColored; //没有被染色的点
  auto & pds = mesh.parallel_data_structure();

  for(auto idx:pds[rank])
    nColored.push_back(idx);

  for(int i=0; i < NE; i++)
    edges.push_back(i);

  int cmax = 0;
  int tnum=0;
  int lnum = nColored.size();
  MPI_Allreduce(&lnum, &tnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  while(tnum > 0)
  {
    cmax++;
     // 只对没有染色的点进行循环
    for(const int& idx: nColored)
    {
      isMin[idx] = true; //初始所有未染色点都是最小点
      randVal[idx] = rand();
    }

    // 通信影像节点上的随机值
    for(const auto & [key, value]:pds)
    {
    }
    

  }

}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream ss;
  ss << argv[1] << "_" << rank << ".vtu";
  std::string fname = ss.str();

  std::cout << "process " << rank << " of " << nprocs << " process." 
    << "The file name is " << fname <<  std::endl;

  PMesh mesh(rank, "node");
  Reader reader(&mesh);
  reader.read(ss.str());

  auto & npid = mesh.node_process_id();
  reader.get_node_data("nid", npid);

  auto & gid = mesh.node_global_id();
  reader.get_node_data("gid", gid);

  mesh.construct_parallel_data_structure();
  std::vector<int> color(); 
  mesh_coloring(mesh, color);

  MPI_Finalize();
  return 0;
}

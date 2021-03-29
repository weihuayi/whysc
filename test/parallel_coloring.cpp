#include <string>
#include <iostream>

#include <mpi.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/VTKMeshReader.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshReader<TriMesh> Reader;


void mesh_coloring(TriMesh & mesh, std::vector<int> & color)
{
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto NN = mesh.number_of_nodes();
  auto NE = mesh.number_of_edges();

  std::vector<int> randVal(NN);
  std::vector<bool> isMin(NN, true);

  std::list<int> edges; //没有被删除的边
  std::list<int> nColored; //没有被染色的点
  auto & pds = mesh.parallel_data_structure();

  for(auto idx:pds[rank])
    nColored.push_back(idx)

  for(int i=0; i < NE; i++)
  {
    edges.push_back(i);
  }

  int cmax = 0;
  int tnum=0;
  int lnum = nColored.size();
  MPI_Allreduce(&lnum, &tnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  while(tnum > 0)
  {
    cmax++;
     // 只对没有染色的点进行循环
    for(const int& index: nColored)
    {
      isMin[index] = true; //初始所有未染色点都是最小点
      randVal[index] = rand();
    }

    // 通信影像节点上的随机值

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
  TriMesh mesh;
  Reader reader(&mesh);
  reader.read(ss.str());

  std::vector<int> nid;
  reader.get_node_data("nid", nid);
  auto & gid = mesh.node_global_id();
  reader.get_node_data("gid", gid);

  auto & pds = mesh.parallel_data_structure();
  for(int j = 0; j < nid.size(); j++)
  {
    pds[nid[j]].push_back(j);
  }

  auto NN = mesh.number_of_nodes();
  std::vector<int> color(NN, 0); 
  mesh_coloring(mesh, color);

  MPI_Finalize();
  return 0;
}

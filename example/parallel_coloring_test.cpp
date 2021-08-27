#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <mpi.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/QuadMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshColoringAlg.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QuadMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, QuadMesh> PMesh;
typedef WHYSC::Mesh::GhostFillingAlg<PMesh> SetGhostAlg;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef WHYSC::Mesh::ParallelMeshColoringAlg<PMesh> PCA;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.get_mesh();

  auto set_ghost_alg = std::make_shared<SetGhostAlg>(mesh, MPI_COMM_WORLD);
  auto color_alg = std::make_shared<PCA>(mesh, set_ghost_alg, MPI_COMM_WORLD);

  auto NN = mesh->number_of_nodes();
  std::vector<int> color(NN);

  color_alg->coloring(true);
  
  MPI_Finalize();
  return 0;
}

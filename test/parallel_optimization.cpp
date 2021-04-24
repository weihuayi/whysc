#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/QuadMesh.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshColoringAlg.h"
#include "mesh/ParallelMeshOptimization.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QuadMesh;
//typedef WHYSC::Mesh::ParallelMesh<GK, QuadMesh> PMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
//typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PMesh;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef WHYSC::Mesh::ParallelMeshColoringAlg<PMesh> PCA;
typedef WHYSC::Mesh::ParallelMeshOptimization<PMesh> PMeshOpt;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::Mesh::EntityOverlap<int> EntityOverlap;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.get_mesh();

  PCA colorAlg(mesh, MPI_COMM_WORLD);

  auto NN = mesh->number_of_nodes();
  std::vector<int> color(NN);

  colorAlg.coloring(color);
  colorAlg.color_test(color);

  PMeshOpt optAlg(mesh, color, MPI_COMM_WORLD);
  for(int i = 0; i < 30; i++)
    optAlg.mesh_optimization();

  std::stringstream ss;
  ss << "opt_"<< mesh->id() << ".vtu";

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();
  writer.write(ss.str());
 
  MPI_Finalize();
  return 0;
}

#include <memory>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/ParallelMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"
#include <iostream>

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;

typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char * argv[])
{
  std::stringstream ss;
  ss << argv[1] << ".vtu";
  auto mesh = std::make_shared<PMesh>();
  Reader reader(mesh);
  reader.read(ss.str());
  auto & gtag = mesh->get_node_int_data()["gtag"];
  auto & gdof = mesh->get_node_int_data()["gdof"];

  reader.get_node_data("gdof", gdof);
  reader.get_node_data("gtag", gtag);

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();
  writer.write("init_tri_mesh.vtu");

  //MF::cube_tetrahedron_mesh(mesh);
  //mesh.uniform_refine(1);
  std::vector<PMesh> submeshes;
  int nparts = std::stoi(argv[2]);

  MF::mesh_node_partition(mesh, nparts, submeshes, "test_tri");
  return 0;
}


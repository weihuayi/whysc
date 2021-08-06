#include <string>
#include <iostream>
#include <vector>

#include "geometry/Geometry_kernel.h"
#include "geometry/CubeWithSpheresModel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/GMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
//typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, HexMesh> PMesh;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::GeometryModel::CubeWithSpheresModel<GK> Model;
typedef WHYSC::Mesh::GMesher<GK, PMesh, Model> GMesher;
typedef WHYSC::Mesh::MeshFactory MF;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;
typedef typename TetMesh::Toplogy Toplogy;

int main()
{
  auto cube = std::make_shared<Model>(2);
  GMesher mesher(cube);
  mesher.mesher3d(0.1, "hex");
  auto mesh = mesher.get_mesh();

  Node n0{1, 2, 3};
  std::cout<< n0[0] << " " << n0[1] << " " << n0[2] <<std::endl;
  //cube->project_to_face(2, n0);
  cube->project_to_face(7, n0);
  std::cout<< n0[0] << " " << n0[1] << " " << n0[2] <<std::endl;
 
  auto & dim = mesh->get_node_int_data()["gdof"];
  auto & tag = mesh->get_node_int_data()["gtag"];

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();

  auto NC = mesh->number_of_cells();
  std::vector<double> z(NC);
  for(int i = 0; i < NC; i++)
  {
    Node n;
    mesh->cell_barycenter(i, n);
    z[i] = n[0]*1000;
  }

  writer.set_point_data(dim, 1, "gdof");
  writer.set_point_data(tag, 1, "gtag");
  writer.set_cell_data(z, 1, "z");
  writer.write("cube.vtu");

  //std::vector<PMesh> submeshes;
  //MF::mesh_node_partition(mesh, 4, submeshes, "test_tet_surface");
  return 0;
}

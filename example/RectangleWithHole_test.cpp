#include <string>
#include <iostream>
#include <vector>

#include "geometry/Geometry_kernel.h"
#include "geometry/RectangleWithHole.h"
#include "mesh/QuadMesh.h"
#include "mesh/TriangleMesh.h"
#include "mesh/GMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/MeshFactory.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, QMesh> PMesh;
typedef WHYSC::Mesh::VTKMeshWriter<QMesh> Writer;
typedef WHYSC::GeometryModel::RectangleWithHole<GK> Model;
typedef WHYSC::Mesh::GMesher<GK, QMesh, Model> GMesher;
typedef WHYSC::Mesh::MeshFactory MF;
typedef WHYSC::Mesh::VTKMeshReader<QMesh> Reader;

/*
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TMesh> PMesh;
typedef WHYSC::Mesh::VTKMeshWriter<TMesh> Writer;
typedef WHYSC::GeometryModel::RectangleWithHole<GK> Model;
typedef WHYSC::Mesh::GMesher<GK, TMesh, Model> GMesher;
typedef WHYSC::Mesh::MeshFactory MF;
typedef WHYSC::Mesh::VTKMeshReader<TMesh> Reader;
*/

int main(int argc, char * argv[])
{
  auto quad = std::make_shared<Model>();
  GMesher mesher(quad);
  std::cout<< argv[1] <<std::endl;
  mesher.mesher2d(0.01, argv[1]);
  auto mesh = mesher.get_mesh();
  auto & nodes = mesh->nodes();

  auto & dim = mesh->get_node_int_data()["gdof"];
  auto & tag = mesh->get_node_int_data()["gtag"];

  int NN = mesh->number_of_nodes();
  for(auto i = 0; i < NN; i++)
  {
    auto & node = nodes[i];
    if(dim[i]==2)
    {
      node[0] += (std::rand()%20 - 10)/3800.0;
    }
  }

  /*
  Node n0{1, 2, 3};
  std::cout<< n0[0] << " " << n0[1] << " " << n0[2] <<std::endl;
  //cube->project_to_face(2, n0);
  cube->project_to_face(7, n0);
  std::cout<< n0[0] << " " << n0[1] << " " << n0[2] <<std::endl;
  */
 

  auto NC = mesh->number_of_cells();
  std::vector<double> z(NC);
  for(int i = 0; i < NC; i++)
  {
    Node n;
    mesh->cell_barycenter(i, n);
    z[i] = n[0]*1000;
  }

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();
  writer.set_point_data(dim, 1, "gdof");
  writer.set_point_data(tag, 1, "gtag");
  writer.set_cell_data(z, 1, "z");
  writer.write("quad_with_quad.vtu");

  std::vector<PMesh> submeshes;
  std::cout<< "hahahah" <<std::endl;
  //MF::mesh_node_partition(mesh, 4, submeshes, "test_tet_surface");
  return 0;
}

#include <string>
#include <iostream>
#include <vector>

#include "geometry/Geometry_kernel.h"
#include "geometry/PolyhedronModel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"
#include "mesh/ParallelMesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
typedef WHYSC::Mesh::VTKMeshWriter<TetMesh> Writer;
typedef WHYSC::GeometryModel::PolyhedronModel<GK, TetMesh> PModel;
typedef WHYSC::Mesh::MeshFactory MF;

int main(int argc, char * argv[])
{
  PModel poly;
  poly.add_point({0.0, 0.0, 0.0}, 1);
  poly.add_point({1.0, 0.0, 0.0}, 2);
  poly.add_point({1.0, 1.0, 0.0}, 3);
  poly.add_point({0.0, 1.0, 0.0}, 4);
  poly.add_point({0.0, 0.0, 1.0}, 5);
  poly.add_point({1.0, 0.0, 1.0}, 6);
  poly.add_point({1.0, 1.0, 1.0}, 7);
  poly.add_point({0.0, 1.0, 1.0}, 8);

  poly.add_line({1, 2}, 1);
  poly.add_line({2, 3}, 2);
  poly.add_line({3, 4}, 3);
  poly.add_line({4, 1}, 4);
  poly.add_line({5, 6}, 5);
  poly.add_line({6, 7}, 6);
  poly.add_line({7, 8}, 7);
  poly.add_line({8, 5}, 8);
  poly.add_line({1, 5}, 9);
  poly.add_line({2, 6}, 10);
  poly.add_line({3, 7}, 11);
  poly.add_line({4, 8}, 12);

  poly.add_face({-1, -4, -3, -2}, 1);
  poly.add_face({5, 6, 7, 8}, 2);
  poly.add_face({1, 10, -5, -9}, 3);
  poly.add_face({3, 12, -7, -11}, 4);
  poly.add_face({2, 11, -6, -10}, 5);
  poly.add_face({4, 9, -8, -12}, 6);

  poly.add_polyhedron({1, 2, 3, 4, 5, 6}, 1);
  std::vector<int> dim;
  std::vector<int> tag;
  poly.mesher(dim, tag);
  auto mesh = poly.get_mesh();

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();

  writer.set_point_data(dim, 1, "dim");
  writer.set_point_data(tag, 1, "tag");
  writer.write("cube.vtu");
  return 0;
}

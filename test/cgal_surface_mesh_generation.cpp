#include <sstream> 

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <metis.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PMesh;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;

FT sphere_function (Point_3 p) 
{
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

FT test_function (Point_3 p)
{
    double x = p.x();
    double y = p.y();
    double z = p.z();
    return pow(x*x+9*y*y/4+z*z-1, 3) - x*x*z*z*z-9*y*y*z*z*z/80;
}


int main() 
{
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  // defining the surface
  Surface_3 surface(test_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 10.)); // bounding sphere
  // Note that "2." above is the *squared* radius of the bounding sphere!
  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.1,  // radius bound
                                                     0.1); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  Surface_mesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

  PMesh mesh(1);
  MF::cgal_surface_mesh_to_triangle_mesh(sm, mesh);

  std::vector<PMesh> submeshes;
  MF::mesh_node_partition(mesh, 4, submeshes, "test_surface");

  /*
  std::vector<std::vector<int>> nids;
  nids.resize(4);
  for(int i = 0; i < 4; i++)
  {
    int NN = submeshes[i].number_of_nodes();
    nids[i].resize(NN);
    for(auto & j : nids[i])
    {
      j = i;
    }
    auto pds = submeshes[i].parallel_data_structure();
    for(int j = 0; j < 4; j++)
    {
      auto it = pds.find(j);
      if(it != pds.end())
      {
        for(auto k : it->second)
        {
          nids[i][k] = j;
        }
      }
    }
    std::stringstream ss;
    ss << "test_surface_" << i << ".vtu";
    Writer writer(&submeshes[i]);
    writer.set_points();
    writer.set_cells();
    writer.set_point_data(nids[i], 1, "nid");
    writer.set_point_data(submeshes[i].node_global_id(), 1, "gid");
    writer.write(ss.str());
  }

  //测试读文件
  std::vector<TriMesh> meshes(4);
  for(int i = 0; i < 4; i++)
  {
    std::stringstream ss;
    ss << "test_surface_" << i << ".vtu";

    Reader reader(&meshes[i]);
    reader.read(ss.str());

    std::vector<int> nid0;
    reader.get_node_data("nid", nid0);

    auto & gid = meshes[i].node_global_id();
    reader.get_node_data("gid", gid);

    auto & pds = meshes[i].parallel_data_structure();
    for(int j = 0; j < nid0.size(); j++)
    {
      if(nid0[j]!=i)
      {
        pds[nid0[j]].push_back(j);
      }
    }
  */
}

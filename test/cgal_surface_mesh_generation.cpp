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
typedef TriMesh::Cell Cell;
typedef TriMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshWriter<TriMesh> Writer;
typedef WHYSC::Mesh::VTKMeshReader<TriMesh> Reader;
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

  TriMesh mesh;
  MF::cgal_surface_mesh_to_triangle_mesh(sm, mesh);

  std::vector<int> nid;
  std::vector<int> cid;
  MF::mesh_node_partition(mesh, 4, nid, cid);

  Writer writer(&mesh);
  writer.set_points();
  writer.set_cells();
  writer.set_point_data(nid, 1, "nid");
  writer.set_cell_data(cid, 1, "cid");
  writer.write("test_surface_0.vtu");

  TriMesh tmesh;
  Reader reader(&tmesh);
  reader.read("test_surface_0.vtu");

  std::vector<int> nid1;
  std::vector<int> cid1;
  reader.get_node_data("nid", nid1);
  reader.get_cell_data("cid", cid1);

  for(auto & i : nid1)
  {
    std::cout<< i <<std::endl;
  }

  for(auto & i : cid1)
  {
    std::cout<< i <<std::endl;
  }

  std::vector<TriMesh> meshs;
  std::vector<std::vector<Node>> nodes;
  std::vector<std::vector<Cell>> cells;
  std::vector<std::map<int, int>> nidxmap;//节点在整体网格中的编号 
  meshs.resize(4);
  cells.resize(4);
  nodes.resize(4);
  nidxmap.resize(4);

  auto & cell = tmesh.cells();
  auto & node = tmesh.cells();

  for(auto & c : cell)
  {
    bool isinmesh[4] = {0};//判断单元在哪个网格中
    for(auto & i : c)
    {
      isinmesh[nid[i]] = true;
    }

    for(int i = 0; i < 4; i++)
    {
      if(isinmesh[i])
      {
        cells[i].push_back(c);
      }
    }
  }

  for(int i = 0; i < 4; i++)
  {
    int num = 0;
    for(auto & c : cells[i])
    {
      for(auto & v : c)
      {
        if(it == nidxmap[i].end())
        {
          nidxmap[i].insert(std::pair<int, int>(v, num));
          nodes[i].push_back(node[v]);
          num++;
        }
      }
    }
  }

  for(int i = 0; i < 4; i++)
  {
    for(auto & c : cells[i])
    {
      for(auto & v : c)
      {
        auto it = nidxmap[i].find(v);
        v = it->second;
      }
    }
  }



  for(int i = 0; i < 4; i++)
  {
    meshs[i].nodes() = nodes[i];
    meshs[i].cells() = cells[i];
    meshs[i].init_top();
  }

  Writer writer0(&meshs[0]);
  writer0.set_points();
  writer0.set_cells();
  writer0.write("test_surface_3.vtu");
  
}

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
#include "mesh/QuadMesh.h"
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
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QuadMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PTMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, QuadMesh> PQMesh;

typedef WHYSC::Mesh::VTKMeshReader<PTMesh> TReader;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::VTKMeshReader<PQMesh> QReader;
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

  PTMesh mesh(1);
  MF::cgal_surface_mesh_to_triangle_mesh(sm, mesh);
  
  //MF::one_triangle_mesh(mesh);
  //mesh.uniform_refine(9);

  std::vector<PTMesh> submeshes;

  //MF::mesh_node_partition(mesh, 4, submeshes, "test_surface");

}

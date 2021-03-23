
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef TetMesh::Cell Cell;
typedef TetMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshWriter<TetMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
typedef CGAL::Sequential_tag Concurrency_tag;

typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

using namespace CGAL::parameters;

FT sphere_function (const Point& p)
{   
    return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; 
}

FT test_function (const Point& p)
{
    double x = p.x();
    double y = p.y();
    double z = p.z();
    return pow(x*x+9*y*y/4+z*z-1, 3) - x*x*z*z*z-9*y*y*z*z*z/80;
}

int main()
{
    Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain(test_function,
                                             K::Sphere_3(CGAL::ORIGIN, 10.));
    // Mesh criteria
    Mesh_criteria criteria(facet_angle=30, facet_size=0.4, facet_distance=0.25,
                         cell_radius_edge_ratio=2, cell_size=0.2);
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    TetMesh mesh;

    MF::c3t3_to_tetmesh(c3t3, mesh);

    Writer writer(&mesh);
    writer.set_points();
    writer.set_cells();
    writer.write("test_surface.vtu");
    return 0;
}


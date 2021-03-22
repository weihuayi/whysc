
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
    double x = p.x();
    double y = p.y();
    double z = p.z();
 
    std::cout<< CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1 <<std::endl;
    std::cout<< "x = " << " " << x <<std::endl; 
    std::cout<< "y = " << " " << y <<std::endl; 
    std::cout<< "y = " << " " << z <<std::endl; 
    return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; 
}

FT test_function (const Point& p)
{
    //double r = p[0]*p[0]*p[0] + p[1]*p[1]*p[1] +p[2]*p[2]*p[2];
    double x = p.x();
    double y = p.y();
    double z = p.z();
    double r = pow(pow(x, 2)+pow(y, 2) + pow(z, 2), 1/2);
    Point pp(CGAL::ORIGIN);
    double x1 = pp.x();
    double y1 = pp.y();
    double z1 = pp.z();
    std::cout<< "r = " << " " << r <<std::endl; 
    std::cout<< "x1 = " << " " << x1 <<std::endl; 
    std::cout<< "y1 = " << " " << y1 <<std::endl; 
    std::cout<< "y1 = " << " " << z1 <<std::endl; 

    Point D(1, 1, 2);
    std::cout<< "test" << " " << CGAL::squared_distance(D, Point(CGAL::ORIGIN)) <<std::endl;
    return r;
}



int main()
{
    Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain(sphere_function,
                                             K::Sphere_3(CGAL::ORIGIN, 2.));
    // Mesh criteria
    Mesh_criteria criteria(facet_angle=30, facet_size=0.8, facet_distance=0.5,
                         cell_radius_edge_ratio=2, cell_size=1);
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


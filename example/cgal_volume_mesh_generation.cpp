
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <memory>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/QuadMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/ParallelMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"
#include <iostream>

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
//typedef WHYSC::Mesh::ParallelMesh<GK, QMesh> PMesh;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;
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

    //PMesh mesh(1);
    //MF::cgal_c3t3_to_tetmesh(c3t3, mesh);

    auto mesh = std::make_shared<PMesh>();
    Reader reader(mesh);
    reader.read("cube.vtu");
    auto & gtag = mesh->get_node_int_data()["gtag"];
    auto & gdof = mesh->get_node_int_data()["gdof"];

    reader.get_node_data("gdof", gdof);
    reader.get_node_data("gtag", gtag);

    //MF::cube_tetrahedron_mesh(mesh);
    //mesh.uniform_refine(1);
    std::vector<PMesh> submeshes;

    MF::mesh_node_partition(mesh, 4, submeshes, "test_tet_surface");

    //std::vector<bool> isBdNode;
    //mesh->is_boundary_node(isBdNode);

    //std::vector<int> fixednode;
    //for(auto i : isBdNode)
    //{
    //  if(i)
    //    fixednode.push_back(1);
    //  else
    //    fixednode.push_back(0);
    //}

    //std::vector<bool> isBdFace;
    //mesh->is_boundary_face(isBdFace);

    //std::vector<int> fixedCell(mesh->number_of_cells());
    //for(int i = 0; i < isBdFace.size(); i++)
    //{
    //  if(isBdFace[i])
    //  {
    //    fixedCell[mesh->face_to_cell(i)[0]] = true;
    //  }
    //}

    //Writer writer(mesh);
    //writer.set_points();
    //writer.set_cells();
    //writer.set_point_data(fixednode, 1, "fixed");
    //writer.set_cell_data(fixedCell, 1, "fixed");
    //writer.write("test_surface.vtu");
    return 0;
}


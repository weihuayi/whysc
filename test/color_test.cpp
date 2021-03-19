
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <metis.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/color.h"

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

int main(int argc, char **argv)
{
   
    Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain(sphere_function,
                                             K::Sphere_3(CGAL::ORIGIN, 2.));
    Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
                         cell_radius_edge_ratio=2, cell_size=0.1);
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    TetMesh mesh;

    MF::c3t3_to_tetmesh(c3t3, mesh);

    //MF::one_tetrahedron_mesh(mesh, "equ");
    //MF::one_tetrahedron_mesh(mesh, "iso");
    //MF::cube_tetrahedron_mesh(mesh);
 
    //
    Toplogy node2node;
    mesh.node_to_node(node2node);

    std::vector<int> &nei = node2node.neighbors();
    std::vector<int> &loc = node2node.locations();

    int NN = mesh.number_of_nodes();
    int col[NN]={0};
    mesh_coloring<TetMesh>(mesh, col);

    //for(int i = 0; i<NN; i++)
    //{
    //    std::cout<< i << "的颜色" << col[i] << std::endl;
    //}

    //for(auto i=mesh.edge_begin(); i!=mesh.edge_end(); i++)
    //{
    //    std::cout<< "edge " << (*i)[0]  << " " << (*i)[1] << std::endl;
    //}

    std::vector<double> nodedata;
    for(int i = 0; i < NN; i++)
    {
        nodedata.push_back((double)col[i]);
    }

    Writer writer(&mesh);
    writer.set_points();
    writer.set_point_data(nodedata, 1, "color");
    writer.set_cells();
    writer.write("color_test.vtu");

    return 0;
}

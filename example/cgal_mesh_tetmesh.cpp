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

int main()
{
    Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain(sphere_function,
                                             K::Sphere_3(CGAL::ORIGIN, 2.));
    // Mesh criteria
    Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
                         cell_radius_edge_ratio=2, cell_size=0.1);
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    TetMesh mesh;

    MF::c3t3_to_tetmesh(c3t3, mesh);

    //Toplogy cell2cell;
    //mesh.cell_to_cell(cell2cell);

    //Toplogy node2node;
    //mesh.node_to_node(node2node);
    //
    Toplogy cell2node;
    mesh.cell_to_node(cell2node);

    idx_t ne = mesh.number_of_cells();
    idx_t nn = mesh.number_of_nodes();
    idx_t ncommon = 3; // 四面体共享三个顶点
    idx_t nparts = 2;
    idx_t objval;
    std::vector<int> epart(ne);
    std::vector<int> npart(nn);

    int r = METIS_PartMeshDual(&ne, &nn, cell2node.locations().data(), cell2node.neighbors().data(),
            NULL, NULL, &ncommon, &nparts, 
            NULL, NULL, &objval, epart.data(), 
            npart.data());

    //int r = METIS_PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
    //                  idx_t *vwgt, idx_t *vsize, idx_t *ncommon, idx_t *nparts, 
    //                 real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, 
    //                 idx_t *npart);

    r = METIS_PartMeshNodal(&ne, &nn, cell2node.locations().data(), cell2node.neighbors().data(),
            NULL, NULL, &nparts, NULL,
            NULL, &objval, epart.data(), npart.data());
    //METIS_API(int) METIS_PartMeshNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
    //                  idx_t *vwgt, idx_t *vsize, idx_t *nparts, real_t *tpwgts, 
    //                 idx_t *options, idx_t *objval, idx_t *epart, idx_t *npart);

    Writer writer(&mesh);
    writer.set_points();
    writer.set_cells();
    writer.set_cell_data(epart, 1, "CID");
    writer.set_point_data(npart, 1, "NID");
    writer.write("test.vtu");
    return 0;
}


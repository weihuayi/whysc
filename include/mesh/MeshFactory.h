#ifndef MeshFactory_h
#define MeshFactory_h

#include <cmath>
#include <string>
#include <unordered_map>
#include <iostream>

namespace WHYSC {
namespace Mesh {

class MeshFactory
{

public:
  template<typename QuadMesh>
  static void one_quad_mesh(QuadMesh & mesh)
  {
    typedef typename QuadMesh::Node Node;
    typedef typename QuadMesh::Cell Cell;
    auto & nodes = mesh.nodes();
    nodes.reserve(4);

    nodes.push_back(Node{0.0, 0.0});
    nodes.push_back(Node{1.0, 0.0});
    nodes.push_back(Node{1.0, 1.0});
    nodes.push_back(Node{0.0, 1.0});

    auto & cells = mesh.cells();
    cells.reserve(1);
    cells.push_back(Cell{0, 1, 2, 3});
    mesh.init_top();
    return;
  }

  template<typename TriMesh>
  static void one_triangle_mesh(TriMesh & mesh, const std::string type="equ")
  {
    typedef typename TriMesh::Node Node;
    typedef typename TriMesh::Cell Cell;
    auto & nodes = mesh.nodes();
    nodes.reserve(3);

    if(type == "equ")
    {
        nodes.push_back(Node{0.0, 0.0});
        nodes.push_back(Node{1.0, 0.0});
        nodes.push_back(Node{0.5, std::sqrt(3.0)/2.0});
    }
    else if(type == "iso")
    {
        nodes.push_back(Node{0.0, 0.0});
        nodes.push_back(Node{1.0, 0.0});
        nodes.push_back(Node{0.0, 1.0});
    }

    auto & cells = mesh.cells();
    cells.reserve(1);
    cells.push_back(Cell{0, 1, 2});
    mesh.init_top();
    return;
  }

  template<typename TriMesh>
  static void square_triangle_mesh(TriMesh & mesh)
  {
    typedef typename TriMesh::Node Node;
    typedef typename TriMesh::Cell Cell;
    mesh.insert(Node{0.0, 0.0});
    mesh.insert(Node{1.0, 0.0});
    mesh.insert(Node{1.0, 1.0});
    mesh.insert(Node{0.0, 1.0});

    mesh.insert(Cell{1, 2, 0});
    mesh.insert(Cell{3, 0, 2});
    mesh.init_top();
    return;
  }

  template<typename TetMesh>
  static void cube_tetrahedron_mesh(TetMesh & mesh)
  {
    typedef typename TetMesh::Node Node;
    typedef typename TetMesh::Cell Cell;
    mesh.insert(Node{0.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 1.0, 1.0});
    mesh.insert(Node{0.0, 1.0, 1.0});

    mesh.insert(Cell{0, 1, 2, 6});
    mesh.insert(Cell{0, 5, 1, 6});
    mesh.insert(Cell{0, 4, 5, 6});
    mesh.insert(Cell{0, 7, 4, 6});
    mesh.insert(Cell{0, 3, 7, 6});
    mesh.insert(Cell{0, 2, 3, 6});
    mesh.init_top();
    return;
  }

  template<typename TetMesh>
  static void one_tetrahedron_mesh(TetMesh & mesh, const std::string type="equ")
  {
    typedef typename TetMesh::Node Node;
    typedef typename TetMesh::Cell Cell;
    auto & nodes = mesh.nodes();
    nodes.reserve(4);

    if(type == "equ")
    {
        nodes.push_back(Node{0.0, 0.0, 0.0});
        nodes.push_back(Node{1.0, 0.0, 0.0});
        nodes.push_back(Node{0.5, std::sqrt(3.0)/2.0, 0.0});
        nodes.push_back(Node{0.5, std::sqrt(3.0)/6.0, std::sqrt(2.0/3.0)}); 
    }
    else if(type == "iso")
    {
        nodes.push_back(Node{0.0, 0.0, 0.0});
        nodes.push_back(Node{1.0, 0.0, 0.0});
        nodes.push_back(Node{0.0, 1.0, 0.0});
        nodes.push_back(Node{0.0, 0.0, 1.0}); 
    }

    auto & cells = mesh.cells();
    cells.reserve(1);
    cells.push_back(Cell{0, 1, 2, 3});
    mesh.init_top();
    return;
  }

  template<typename C3T3,  typename TetMesh>
  static void c3t3_to_tetmesh(C3T3 & c3t3, TetMesh & mesh)
  {
    typedef typename C3T3::Triangulation Tr;
    typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
    typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Point Point; //can be weighted or not

    typedef typename TetMesh::Node Node;

    const Tr& tr = c3t3.triangulation();

    auto & nodes = mesh.nodes();
    nodes.resize(tr.number_of_vertices());

    std::unordered_map<Vertex_handle, int> V;
    int i = 0;
    for(auto vit = tr.finite_vertices_begin();
            vit != tr.finite_vertices_end();
            ++vit, ++i)
    {
        V[vit] = i;
        auto p = tr.point(vit);
        nodes[i][0] = p.x();
        nodes[i][1] = p.y();
        nodes[i][2] = p.z();
    }

    auto & cells = mesh.cells();
    cells.resize(c3t3.number_of_cells_in_complex());
    i = 0;
    for( auto cit = c3t3.cells_in_complex_begin();
       cit != c3t3.cells_in_complex_end();
       ++cit, ++i)
    {
        for(int j=0; j<4; j++)
          cells[i][j] = V[cit->vertex(j)];
    }
    mesh.init_top();
  }

  template<typename C2T3,  typename TriMesh>
  static void C2T3_to_triangle_mesh(C2T3 & c2t3, TriMesh & mesh)
  {
    typedef typename C2T3::Triangulation Tr;

    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Point Point; //can be weighted or not

    typedef typename TriMesh::Node Node;

    const Tr& tr = c2t3.triangulation();

    auto & nodes = mesh.nodes();
    nodes.resize(tr.number_of_vertices());

    std::unordered_map<Vertex_handle, int> V;
    int i = 0;
    for(auto vit = tr.finite_vertices_begin();
            vit != tr.finite_vertices_end();
            ++vit, ++i)
    {
        std::cout<< i <<std::endl;
        V[vit] = i;
        auto p = tr.point(vit);
        nodes[i][0] = p.x();
        nodes[i][1] = p.y();
        nodes[i][2] = p.z();
    }

    auto & cells = mesh.cells();
    cells.resize(c2t3.number_of_facets());
    i = 0;
    for( auto fit = c2t3.facets_begin();
       fit != c2t3.facets_end();
       ++fit, ++i)
    {
      for(int j=0; j<3; j++)
        cells[i][j] = V[fit->first->vertex(j)];
    }
    mesh.init_top();

  }
  template<typename SMesh,  typename TriMesh>
  static void Surface_mesh_to_triangle_mesh(SMesh & sm, TriMesh & mesh)
  {
    typedef typename SMesh::Vertex_index Vindex;
    auto & nodes = mesh.nodes();
    nodes.resize(sm.number_of_vertices());

    std::unordered_map<Vindex, int> V;
    int i = 0;
    auto vr = sm.vertices();
    for(auto vit = vr.begin(); vit != vr.end(); ++vit, ++i)
    {
        V[*vit] = i;
        auto p = sm.point(*vit);
        nodes[i][0] = p.x();
        nodes[i][1] = p.y();
        nodes[i][2] = p.z();
    }

    auto & cells = mesh.cells();
    cells.resize(sm.number_of_faces());
    i = 0;
    auto fr = sm.faces();
    for( auto fit = fr.begin(); fit != fr.end(); ++fit, ++i)
    {
      auto h0 = sm.halfedge(*fit);
      auto h1 = sm.next(h0);
      auto h2 = sm.next(h1);

      cells[i][0] = V[sm.target(h0)];
      cells[i][1] = V[sm.target(h1)];
      cells[i][2] = V[sm.target(h2)];
    }
    mesh.init_top();
  }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshFactory_h

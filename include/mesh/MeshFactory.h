#ifndef MeshFactory_h
#define MeshFactory_h

#include <cmath>
#include <string>
#include <unordered_map>
#include <iostream>
#include <set>

#include <metis.h>

#include "VTKMeshWriter.h"

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
  static void cgal_c3t3_to_tetmesh(C3T3 & c3t3, TetMesh & mesh)
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

  template<typename SMesh,  typename TriMesh>
  static void cgal_surface_mesh_to_triangle_mesh(SMesh & sm, TriMesh & mesh)
  {
    typedef typename SMesh::Vertex_index VIndex;
    auto & nodes = mesh.nodes();
    nodes.resize(sm.number_of_vertices());

    //std::unordered_map<VIndex, int> V;
    for(auto & v : sm.vertices())
    {
        auto i = v.idx();
        auto p = sm.point(v);
        nodes[i][0] = p.x();
        nodes[i][1] = p.y();
        nodes[i][2] = p.z();
    }

    auto & cells = mesh.cells();
    cells.resize(sm.number_of_faces());
    auto fr = sm.faces();

    for( auto & f : sm.faces())
    {
      auto i = f.idx();
      auto h0 = sm.halfedge(f);
      auto h1 = sm.next(h0);
      auto h2 = sm.next(h1);

      cells[i][0] = sm.target(h0).idx();
      cells[i][1] = sm.target(h1).idx();
      cells[i][2] = sm.target(h2).idx();
    }

    mesh.init_top();
  }

  template<typename Mesh>
  static void mesh_node_partition(Mesh & mesh, idx_t nparts, std::vector<int> & nid, std::vector<int> & cid)
  {
    typedef typename Mesh::Toplogy Toplogy;

    Toplogy cell2node;
    mesh.cell_to_node(cell2node);

    idx_t ne = mesh.number_of_cells();
    idx_t nn = mesh.number_of_nodes();
    cid.resize(ne);
    nid.resize(nn);

    idx_t objval;
    auto r = METIS_PartMeshNodal(&ne, &nn, cell2node.locations().data(), cell2node.neighbors().data(),
            NULL, NULL, &nparts, NULL,
            NULL, &objval, cid.data(), nid.data());
  }

  template<typename Mesh>
  static void mesh_node_partition(Mesh & mesh, idx_t nparts, std::vector<Mesh> & submeshes, 
      std::string fname="")
  {
    typedef typename Mesh::Toplogy Toplogy;
    typedef VTKMeshWriter<Mesh> Writer;

    Toplogy cell2node;
    mesh.cell_to_node(cell2node);

    idx_t ne = mesh.number_of_cells();
    idx_t nn = mesh.number_of_nodes();

    std::vector<int> cid(ne);
    std::vector<int> nid(nn);

    idx_t objval;
    auto r = METIS_PartMeshNodal(&ne, &nn, cell2node.locations().data(), cell2node.neighbors().data(),
            NULL, NULL, &nparts, NULL,
            NULL, &objval, cid.data(), nid.data());
    submeshes.resize(nparts);


    auto & cells = mesh.cells();
    auto & nodes = mesh.nodes();

    //将 cell 分到每个 mesh
    for(auto & cell : cells)
    {
      std::set<int> idx; // 存放 cell 都在哪个网格中
      for(auto i: cell)
      {
        idx.insert(nid[i]); // cell 在其节点所在网格
      }

      for(auto i: idx)
      {
        submeshes[i].cells().push_back(cell); // 将 cell 放到其所在网格
      }
    }

    std::vector<int> nums(nparts, 0);
    for(int gid = 0; gid < nn; gid++)
    {
      submeshes[nid[gid]].nodes().push_back(nodes[gid]); // 做映射
      submeshes[nid[gid]].node_global_id().push_back(gid);

      submeshes[nid[gid]].node_global_to_local_id().insert(std::pair<int, int>(gid, nums[nid[gid]]));
      nums[nid[gid]] += 1;
    }

    for(int i = 0; i < nparts; i++)
    {
      auto & pds = submeshes[i].parallel_data_structure();
      auto & ng2l = submeshes[i].node_global_to_local_id();

      for(auto & cell: submeshes[i].cells())
      {
        for(auto & v : cell)
        {
          auto it = ng2l.find(v);
          if(it == ng2l.end())
          {
            submeshes[i].nodes().push_back(nodes[v]);
            submeshes[i].node_global_id().push_back(v);
            ng2l.insert(std::pair<int, int>(v, nums[i]));
            pds[nid[v]].insert(nums[i]);
            v = nums[i];
            nums[i] += 1;
          }
          else
          {
            v = it->second; 
          }
        }
      }
    }

    if(!fname.empty())
    {
      std::vector<std::vector<int>> nids;
      nids.resize(nparts);
      for(int i = 0; i < nparts; i++)
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
        ss << fname <<"_"<< i << ".vtu";
        Writer writer(&submeshes[i]);
        writer.set_points();
        writer.set_cells();
        writer.set_point_data(nids[i], 1, "nid");
        writer.set_point_data(submeshes[i].node_global_id(), 1, "gid");
        writer.write(ss.str());
      }
    }
  }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshFactory_h

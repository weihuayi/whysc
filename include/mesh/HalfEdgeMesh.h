#ifndef HalfEdgeMesh_h
#define HalfEdgeMesh_h

#include <vector>
#include <array>
#include <map>

#include "MeshToplogy.h"
#include "thirdparty/json.hpp"

using json = nlohmann::json;

namespace WHYSC {
namespace Mesh {

/*
 * 
 *
 * Notes
 * -----
 *  四面体网格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR>
class HalfEdgeMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;
  typedef typename GK::Int I;
  typedef typename GK::Float F;

  // 实体类型
  typedef typename std::array<I, 2> Edge;
  typedef typename std::array<I, 5> HalfEdge;
  typedef typename std::array<I, 4> Face;
  typedef typename std::array<I, 8> Cell;

  typedef std::vector<Node> NodeArray;
  typedef std::vector<Edge> EdgeArray;
  typedef std::vector<HalfEdge> HalfEdgeArray;
  typedef std::vector<Face> FaceArray;
  typedef std::vector<Cell> CellArray;

  // 迭代子类型
  typedef typename NodeArray::iterator NodeIterator;
  typedef typename NodeArray::const_iterator ConstNodeIterator;

  typedef typename EdgeArray::iterator EdgeIterator;
  typedef typename EdgeArray::const_iterator ConstEdgeIterator;

  typedef typename FaceArray::iterator FaceIterator;
  typedef typename FaceArray::const_iterator ConstFaceIterator;

  typedef typename CellArray::iterator CellIterator;
  typedef typename CellArray::const_iterator ConstCellIterator;

public:

  HalfEdgeMesh()
  {
  }

  void insert(const Node & node)
  {
    m_node.push_back(node);
  }

  void insert(const HalfEdgeMesh & halfedge)
  {
    m_halfedge.push_back(halfedge);
  }

  I number_of_nodes()
  {
    return m_node.size();
  }

  I number_of_halfedges()
  {
    return m_halfedge.size();
  }

  static I geo_dimension()//TODO
  {
      return 3;
  }

  static I top_dimension()//TODO
  {
      return 3;
  }

  I vtk_cell_type(I TD=3)//TODO
  {
      if(TD == 3)
          return 72; // VTK_TETRA = 10
      else if(TD == 2)
          return 70; // VTK_TRIANGLE = 5
      else if(TD == 1)
          return 3; // VTK_LINE = 1
      else
          return 68; // VTK_EMPLTY_CELL = 0
  }

  /*
   * 
   * Notes
   * -----
   *  TODO: 考虑如何在原来拓扑的基础上重建拓扑信息
   */
  void update_top()
  {
      return;
  }

  template<typename Mesh>
  void from_mesh(Mesh & mesh)
  {
    m_node = mesh.nodes();

    auto NN = mesh.number_of_nodes();
    auto NE = mesh.number_of_edges();
    auto & edge = mesh.edges();

    m_halfedge.resize(NE*2);
    for(int i = 0; i < NE; i++)
    {
      auto & e2c = mesh.edge_to_cell(i);

      m_halfedge[2*i][0] = edge[i][0];
      m_halfedge[2*i][1] = e2c[0]+1;//所有边在相邻的第 0 个单元中逆时针
      m_halfedge[2*i][4] = 2*i+1;
      m_halfedge[2*i+1][0] = edge[i][1];
      m_halfedge[2*i+1][4] = 2*i;
      if(e2c[0]==e2c[1])
        m_halfedge[2*i+1][1] = 0;//先把所有的边界半边的单元设为 0, 后面再改变.
      else
        m_halfedge[2*i+1][1] = e2c[1]+1;
    }

    std::map<I, I> idxmap;
    for(int i = 0; i < NE*2; i++)
    {
      auto & h0 = m_halfedge[i];
      auto & h1 = m_halfedge[h0[4]];
      int s = h1[0]+(h0[1]+NN)*(h0[1]+NN+1)/2;
      idxmap[s] = i;
    }

    for(int i = 0; i < NE*2; i++)
    {
      auto & h = m_halfedge[i];
      int s = h[0]+(h[1]+NN)*(h[1]+NN+1)/2;
      auto it = idxmap.find(s);
      m_halfedge[i][2] = it->second;
    }

    for(int i = 0; i < NE*2; i++)
    {
      auto & h = m_halfedge[i];
      m_halfedge[h[2]][3] = i;
    }
  }

  NodeIterator node_begin()
  {
      return m_node.begin();
  }

  NodeIterator node_end()
  {
      return m_node.end();
  }

  NodeArray & nodes()
  {
      return m_node;
  }

  HalfEdgeArray & halfedges()
  {
      return m_halfedge;
  }

  Node & node(const I i)
  {
      return m_node[i];
  }

  HalfEdge & halfedge(const I i)
  {
    return m_halfedge[i];
  }

  json & data()
  {
    return m_data;
  }
  // 实体测度 

  // 实体重心
  void is_boundary_edge(std::vector<bool> & isBdFace)
  {
      auto NHE = number_of_halfedges();
      isBdFace.resize(NHE);
      for(int i = 0; i < NHE; i++)
      {
        if(m_halfedge[i][1] < m_cellstart)
          isBdFace[i] = false;
        else
          isBdFace[i] = true;
      }
  }

  void is_boundary_node(std::vector<bool> & isBdNode)
  {
      auto NN = number_of_nodes();
      auto NHE = number_of_halfedges();
      isBdNode.resize(NN);
      for(int i = 0; i < NHE; i++)
      {
        if(m_halfedge[i][1] < m_cellstart)
          isBdNode[m_halfedge[i][0]] = false;
        else
          isBdNode[m_halfedge[i][0]] = true;
      }
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      print_entities(m_node);

      std::cout << "HalfEdges:" << std::endl;
      print_entities(m_halfedge);
  }

  template<typename Entities>
  void print_entities(Entities & entities)
  {
      auto N = entities.size();
      for(I i = 0; i < N; i++)
      {
          auto & e = entities[i];
          auto n = e.size();
          std::cout << i << ":";
          for(I j = 0; j < n; j++)
          {
              std::cout << " " << e[j];
          }
          std::cout << std::endl;
      }
  }

private:
    int m_cellstart;
    NodeArray m_node;
    HalfEdgeArray m_halfedge;
    std::vector<I> m_hedge;
    std::vector<I> m_hcell;
    json m_data;
};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of HalfEdgeMesh_h

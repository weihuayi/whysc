#ifndef TriangleMesh_h
#define TriangleMesh_h

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
 *  三角形网格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体, 这里 face 实体即为 edge 实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR>
class TriangleMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;

  typedef typename GK::Int I;
  typedef typename GK::Float F;

  typedef typename std::array<I, 3> Cell;
  typedef typename std::array<I, 2> Edge;
  typedef Edge Face;

  typedef typename std::array<I, 4> Edge2cell;
  typedef typename std::array<I, 4> Face2cell;
  typedef typename std::array<I, 3> Cell2edge;
  typedef typename std::array<I, 3> Cell2face;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I> Toplogy;

  typedef std::vector<Node> NodeArray;
  typedef std::vector<Edge> EdgeArray;
  typedef std::vector<Face> FaceArray;
  typedef std::vector<Cell> CellArray;

  typedef typename NodeArray::iterator NodeIterator;
  typedef typename NodeArray::const_iterator ConstNodeIterator;

  typedef typename EdgeArray::iterator EdgeIterator;
  typedef typename EdgeArray::const_iterator ConstEdgeIterator;

  typedef typename FaceArray::iterator FaceIterator;
  typedef typename FaceArray::const_iterator ConstFaceIterator;

  typedef typename CellArray::iterator CellIterator;
  typedef typename CellArray::const_iterator ConstCellIterator;


public:
  TriangleMesh()
  {
      m_holes = 1; // 默认一个外部无界区域的洞
      m_genus = 0;
  }

  void insert(const Node & node)
  {
      m_node.push_back(node);
  }

  void insert(const Cell & cell)
  {
      m_cell.push_back(cell);
  }

  int & number_of_holes()
  {
      return m_holes;
  }

  int & number_of_genus()
  {
      return m_genus;
  }

  static int number_of_nodes_of_each_cell()
  {
      return 3;
  }

  static int number_of_vertices_of_each_cell()
  {
      return 3;
  }

  static int number_of_edges_of_each_cell()
  {
      return 3;
  }

  static int number_of_faces_of_each_cell()
  {
      return 3;
  }

  I number_of_nodes()
  {
      return m_node.size();
  }

  I number_of_cells()
  {
      return m_cell.size();
  }

  I number_of_edges()
  {
      return m_edge.size();
  }

  I number_of_faces()
  {
      return m_edge.size();
  }

  I vtk_cell_type(int TD=2)
  {
      if(TD == 2)
          return 5; // VTK_TRIANGLE
      else if(TD == 1)
          return 3; // VTK_LINE
      else
        return 0;// VTK_EMPLTY_CELL = 0
  }

  int* vtk_cell_index()
  {
    return m_vtkidx;
  }

  int geo_dimension()
  {
      return Node::dimension();
  }

  static int top_dimension()
  {
      return 2;
  }

  /*
   * 
   * Notes
   * -----
   *  TODO: 考虑如何在原来拓扑的基础上重建拓扑信息
   *        原来的边每个变成 2 条, 每个单元内部增加 3 条边
   *        每个单元变成 4 个单元, 这样会提高程序的效率吗?
   */
  void update_top()
  {
      return;
  }

  void init_top()
  {
      auto NN = number_of_nodes();
      auto NC = number_of_cells();
      m_cell2edge.resize(m_cell.size());
      m_edge2cell.clear();
      // 在知道网格代表曲面的洞和亏格的个数后, 可以准确计算边的个数
      m_edge2cell.reserve(NN + NC + m_holes + 2*m_genus - 2);
      std::map<I, I> idxmap;

      I NE = 0;
      // 偏历所有单元
      for(I i = 0; i < m_cell.size(); i++)
      {
          for(I j = 0; j < 3; j++)
          {
             auto s = local_edge_index(i, j);
             auto it = idxmap.find(s);
             if(it == idxmap.end())
             {
                m_cell2edge[i][j] = NE;
                idxmap.insert(std::pair<I, I>(s, NE));
                m_edge2cell.push_back(Edge2cell{i, i, j, j});
                NE++;
             }
             else
             {
                m_cell2edge[i][j] = it->second;
                m_edge2cell[it->second][1] = i;
                m_edge2cell[it->second][3] = j;
             }
          }
      }

      m_edge.resize(NE);
      for(I i = 0; i < NE; i++)
      {
          auto & c = m_cell[m_edge2cell[i][0]];
          auto j = m_edge2cell[i][2];
          m_edge[i][0] = c[m_localedge[j][0]];
          m_edge[i][1] = c[m_localedge[j][1]];
      }
  }

  void is_boundary_edge(std::vector<bool> & isBdEdge)
  {
      auto NE = number_of_edges();
      isBdEdge.resize(NE);
      for(int i = 0; i < NE; i++)
      {
          isBdEdge[i] = false;
          if(edge_to_cell(i)[0]==edge_to_cell(i)[1])
              isBdEdge[i] = true;
      }
  }

  void is_boundary_node(std::vector<bool> & isBdNode)
  {
      auto NN = number_of_nodes();
      auto NE = number_of_edges();
      isBdNode.resize(NN);
      for(int i = 0; i < NE; i++)
      {
          if(edge_to_cell(i)[0]==edge_to_cell(i)[1])
          {
              isBdNode[edge(i)[0]] = true;
              isBdNode[edge(i)[1]] = true;
          }
      }
  }


  void cell_to_cell(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NE = number_of_edges();
      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      for(I i=0; i < NE; i++)
      {
          loc[m_edge2cell[i][0]+1] += 1;
          if(m_edge2cell[i][0] != m_edge2cell[i][1])
          {
              loc[m_edge2cell[i][1]+1] += 1;
          }
      }
      for(I i=0; i < NC; i++)
      {
          loc[i+1] += loc[i];
      }

      auto & nei = top.neighbors();
      nei.resize(loc[NC]);
      std::vector<I> start = loc;
      for(I i = 0; i < NE; i++)
      {
          nei[start[m_edge2cell[i][0]]++] = m_edge2cell[i][1];
          if(m_edge2cell[i][0] != m_edge2cell[i][1])
          {
              nei[start[m_edge2cell[i][1]]++] = m_edge2cell[i][0];
          }
      }
  }

  void node_to_node(Toplogy & top)
  {
    auto NN = number_of_nodes();
    auto NE = number_of_edges();
    auto & loc = top.locations();
    loc.resize(NN+1, 0);
    for(I i=0; i < NE; i++)
    {
      loc[m_edge[i][0]+1] += 1;
      loc[m_edge[i][1]+1] += 1;
    }
    for(I i=0; i < NN; i++)
    {
      loc[i+1] += loc[i];
    }

    auto & nei = top.neighbors();
    nei.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NE; i++)
    {
      nei[start[m_edge[i][0]]++] = m_edge[i][1];
      nei[start[m_edge[i][1]]++] = m_edge[i][0];
    }
  }

  void node_to_cell(Toplogy & top)
  {
    auto NN = number_of_nodes();
    auto NC = number_of_cells();

    auto & loc = top.locations();
    loc.resize(NN+1, 0);
    for(I i=0; i < NC; i++)
    {
      for(auto v : m_cell[i])
        loc[v+1] += 1;
    }
    for(I i=0; i < NN; i++)
    {
      loc[i+1] += loc[i];
    }

    auto & nei = top.neighbors();
    nei.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NC; i++)
    {
      for(auto v : m_cell[i])
        nei[start[v]++] = i;
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

  EdgeIterator edge_begin()
  {
    return m_edge.begin();
  }

  EdgeIterator edge_end()
  {
    return m_edge.end();
  }

  FaceIterator face_begin()
  {
    return m_edge.begin();
  }

  FaceIterator face_end()
  {
    return m_edge.end();
  }

  CellIterator cell_begin()
  {
    return m_cell.begin();
  }

  CellIterator cell_end()
  {
    return m_cell.end();
  }


  std::vector<Node> & nodes()
  {
    return m_node;
  }

  std::vector<Cell> & cells()
  {
    return m_cell;
  }

  std::vector<Edge> & edges()
  {
    return m_edge;
  }

  std::vector<Face> & faces()
  {
    return m_edge;
  }

  Node & node(const I i)
  {
    return m_node[i];
  }

  Cell & cell(const I i)
  {
    return m_cell[i];
  }

  Edge & edge(const I i)
  {
    return m_edge[i];
  }

  Edge2cell & edge_to_cell(const I i)
  {
    return m_edge2cell[i];
  }

  Cell2edge & cell_to_edge(const I i)
  {
    return m_cell2edge[i];
  }

  json & data()
  {
    return m_data;
  }

  F cell_measure(const I i)
  {//TODO: 需要考虑 inline 函数吗?
    auto & c = m_cell[i];
    auto v1 = m_node[c[1]] - m_node[c[0]];
    auto v2 = m_node[c[2]] - m_node[c[0]];
    return 0.5*cross(v1, v2);
  }

  void cell_measure(std::vector<F> & measure)
  {
    auto NC = number_of_cells();
    measure.resize(NC);
    for(I i = 0; i < NC; i++)
        measure[i] = cell_measure(i);
  }

  Node edge_barycenter(const I i)
  {
    auto & e = m_edge[i];
    F x = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
    F y = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
    return Node(x, y);
  }

  void edge_barycenter(const I i, Node & node)
  {
    auto & e = m_edge[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[e[0]][i] + m_node[e[1]][i])/2.0;
  }

  void cell_barycenter(const I i, Node & node)
  {
    auto & c = m_cell[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[c[0]][i] + m_node[c[1]][i] + m_node[c[2]][i])/3.0;
  }

  F edge_measure(const I i)
  {
    auto & e = m_edge[i];
    auto v = m_node[e[1]] - m_node[e[0]];
    return std::sqrt(v.squared_length());
  }

  void edge_measure(std::vector<F> & measure)
  {
    auto NE = number_of_edges();
    measure.resize(NE);
    for(I i = 0; i < NE; i++)
      measure[i] = edge_measure(i);
  }

  F face_measure(const I i)
  {
    auto & e = m_edge[i];
    auto v = m_node[e[1]] - m_node[e[0]];
    return std::sqrt(v.squared_length());
  }

  void face_measure(std::vector<F> & measure)
  {
    auto NE = number_of_edges();
    measure.resize(NE);
    for(I i = 0; i < NE; i++)
      measure[i] = edge_measure(i);
  }

  void uniform_refine(const int n=1)
  {
    for(I i=0; i < n; i++)
    {
      auto NN = number_of_nodes();
      auto NE = number_of_edges();
      m_node.resize(NN + NE);
      for(I j = 0; j < NE; j++)
      {
         edge_barycenter(j, m_node[NN+j]); 
      }
      auto NC = number_of_cells();
      m_cell.resize(4*NC);
      for(I j = 0; j < NC; j++)
      { //TODO: 考虑不同的排列顺序是否影响程序的效率
        auto c = m_cell[j]; 
        m_cell[j][0] = c[0];
        m_cell[j][1] = m_cell2edge[j][2] + NN;
        m_cell[j][2] = m_cell2edge[j][1] + NN;

        m_cell[NC + j][0] = c[1];
        m_cell[NC + j][1] = m_cell2edge[j][0] + NN;
        m_cell[NC + j][2] = m_cell2edge[j][2] + NN;

        m_cell[2*NC + j][0] = c[2];
        m_cell[2*NC + j][1] = m_cell2edge[j][1] + NN;
        m_cell[2*NC + j][2] = m_cell2edge[j][0] + NN;

        m_cell[3*NC + j][0] = m_cell2edge[j][0] + NN; 
        m_cell[3*NC + j][1] = m_cell2edge[j][1] + NN;
        m_cell[3*NC + j][2] = m_cell2edge[j][2] + NN;
      }
      m_edge.clear();
      m_cell2edge.clear();
      m_edge2cell.clear();
      init_top(); 
    }
  }

  void cell_to_node(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NN = number_of_nodes();
      auto nn = number_of_vertices_of_each_cell();

      top.init(2, 0, NC, NN);

      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      auto & nei = top.neighbors();
      nei.resize(NC*nn);

      for(I i = 0; i < NC; i++)
      {
          loc[i+1] = loc[i] + nn;
          nei[nn*i] = m_cell[i][0];
          nei[nn*i + 1] = m_cell[i][1];
          nei[nn*i + 2] = m_cell[i][2];
      }
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      print_entities(m_node);

      std::cout << "Edges:" << std::endl;
      print_entities(m_edge);

      std::cout << "Cells:" << std::endl;
      print_entities(m_cell);

      std::cout << "Edge2cell:" << std::endl;
      print_entities(m_edge2cell);

      std::cout << "Cell2edge:" << std::endl;
      print_entities(m_cell2edge);
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
  /*
   *
   * Notes
   *  计算第 i 个单元的第 j 条边全局唯一的一个整数索引
   */
  I local_edge_index(I i, I j)
  {
      I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
      std::sort(e, e+2);
      return  e[0] + e[1]*(e[1]+1)/2;
  }

private:
  static int m_localedge[3][2];
  static int m_localface[3][2];
  static int m_vtkidx[3];
  int m_holes; // 网格中洞的个数
  int m_genus; // 网格表示曲面的亏格数
  std::vector<Node> m_node;
  std::vector<Edge> m_edge;
  std::vector<Cell> m_cell; 
  std::vector<Edge2cell> m_edge2cell;
  std::vector<Cell2edge> m_cell2edge;
  json m_data;
};


template<typename GK, typename NODE, typename VECTOR>
int TriangleMesh<GK, NODE, VECTOR>::m_localedge[3][2] = {
    {1, 2}, {2, 0}, {0, 1}
};

template<typename GK, typename NODE, typename VECTOR>
int TriangleMesh<GK, NODE, VECTOR>::m_localface[3][2] = {
    {1, 2}, {2, 0}, {0, 1}
};

template<typename GK, typename NODE, typename VECTOR>
int TriangleMesh<GK, NODE, VECTOR>::m_vtkidx[3] = {0, 1, 2};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of TriangleMesh_h

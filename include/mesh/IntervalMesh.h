#ifndef IntervalMesh_h
#define IntervalMesh_h

#include <vector>
#include <array>
#include <map>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "MeshToplogy.h"

namespace WHYSC {
namespace Mesh {

/*
 * 
 *
 * Notes
 * -----
 * 区间单元网格格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR>
class IntervalMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;
  typedef typename GK::Int I;
  typedef typename GK::Float F;

  // 实体类型
  typedef typename std::array<I, 2> Edge;
  typedef typename std::array<I, 2> Face;
  typedef typename std::array<I, 2> Cell;

  typedef typename std::array<I, 2> Node2cell;

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

  static int m_vtkidx[2];

public:
  IntervalMesh()
  {
  }

  void insert(const Node & node)
  {
    m_node.push_back(node);
  }

  void insert(const Cell & cell)
  {
    m_cell.push_back(cell);
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
    return m_cell.size();
  }

  I number_of_faces()
  {
    return m_cell.size();
  }

  static I number_of_nodes_of_each_cell()
  {
      return 2;
  }

  static I number_of_vertices_of_each_cell()
  {
      return 2;
  }

  static I geo_dimension()
  {
      return Node::dimension();
  }

  static I top_dimension()
  {
      return 1;
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
    return m_cell.begin();
  }

  EdgeIterator edge_end()
  {
    return m_cell.end();
  }

  FaceIterator face_begin()
  {
    return m_cell.begin();
  }

  FaceIterator face_end()
  {
    return m_cell.end();
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

  std::vector<Node> & edges()
  {
      return m_cell;
  }

  std::vector<Cell> & faces()
  {
      return m_cell;
  }

  std::vector<Cell> & cells()
  {
      return m_cell;
  }

  Node & node(const I i)
  {
      return m_node[i];
  }

  Edge & edge(const I i)
  {
      return m_cell[i];
  }

  Face & face(const I i)
  {
      return m_cell[i];
  }

  Cell & cell(const I i)
  {
      return m_cell[i];
  }

  I vtk_cell_type(int TD=1)
  {
    return 3; // VTK_LINE
  }

  int* vtk_write_cell_index()
  {
    return m_vtkidx;
  }

  int* vtk_read_cell_index()
  {
    return m_vtkidx;
  }

  int* abaqus_cell_index()
  {
    return m_vtkidx;
  }

  /*
   * 初始化节点和单元的邻接关系信息
   */
  void init_top()
  {
    auto NN = number_of_nodes();
    auto NC = number_of_cells();
    m_node2cell.resize(NN);
    std::vector<int> start(NN, 0);

    // 遍历所有单元
    for(I i = 0; i < NC; i++)
    {
      for(auto v : m_cell[i])
      {
        m_node2cell[v][start[v]] = i;
        start[v]++;
      }
    }
  }

  void uniform_refine(int n = 1)
  {
    for(int it = 0; it < n; it++)
    {
      int NC = number_of_cells();
      int NN = number_of_nodes();
      m_node.resize(NN+NC);
      for(int i = 0; i < NC; i++)
      {
        for(int j = 0; j < geo_dimension(); j++)
        {
          m_node[NN+i][j] = (m_node[m_cell[i][0]][j] + m_node[m_cell[i][1]][j])/2;
        }
      }
      m_cell.resize(NC*2);
      for(int i = 0; i < NC; i++)
      {
        m_cell[NC+i][0] = NN+i;
        m_cell[NC+i][1] = m_cell[i][1];
        m_cell[i][1] = NN+i;
      }
      init_top();
    }
  }

  void print()
  {
    std::cout << "Nodes:" << std::endl;
    print_entities(m_node);

    std::cout << "Cells:" << std::endl;
    print_entities(m_cell);
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
     * -----
     *  第 i 个 cell 的第 j 个 node 的全局唯一整数编号
     */
    I local_node_index(I i, I j)
    {
        return  m_cell[i][j];
    }

private:
    std::vector<Node> m_node;
    std::vector<Cell> m_cell; 

    std::vector<Node2cell> m_node2cell;
};

template<typename GK, typename Node, typename Vector>
int IntervalMesh<GK, Node, Vector>::m_vtkidx[2] = {0, 1};
} // end of namespace Mesh 

} // end of namespace TOPT

#endif // end of IntervalMesh_h

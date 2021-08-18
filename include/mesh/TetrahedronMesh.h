#ifndef TetrahedronMesh_h
#define TetrahedronMesh_h

#include <vector>
#include <array>
#include <map>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "MeshToplogy.h"
#include "NodeData.h"

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
class TetrahedronMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;
  typedef typename GK::Int I;
  typedef typename GK::Float F;

  // 实体类型
  typedef typename std::array<I, 2> Edge;
  typedef typename std::array<I, 3> Face;
  typedef typename std::array<I, 4> Cell;

  // 规则化的实体关系类型，这里的规则化是指一种实体有固定个数的另一种邻接实体
  // 每个单元有 4 个三角形面
  // 每个单元有 6 条边
  typedef typename std::array<I, 4> Face2cell; // 如果边界面则左右单元的编号一样
  typedef typename std::array<I, 4> Cell2face;
  typedef typename std::array<I, 4> Cell2cell;
  typedef typename std::array<I, 6> Cell2edge;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I, std::vector<I> > Toplogy;

  // 迭代子类型
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

  static int m_localedge[6][2];
  static int m_localface[4][3];
  static int m_localface2edge[4][3];
  static int m_refine[3][6];
  static int m_index[12][4];
  static int m_vtkidx[4];
  static int m_num[4][4];

public:

  TetrahedronMesh()
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
    return m_edge.size();
  }

  I number_of_faces()
  {
    return m_face.size();
  }

  static I number_of_nodes_of_each_cell()
  {
      return 4;
  }

  static I number_of_vertices_of_each_cell()
  {
      return 4;
  }

  static I geo_dimension()
  {
      return 3;
  }

  static I top_dimension()
  {
      return 3;
  }

  int* vtk_read_cell_index()
  {
    return m_vtkidx;
  }

  int* vtk_write_cell_index()
  {
    return m_vtkidx;
  }


  I vtk_cell_type(I TD=3)
  {
      if(TD == 3)
          return 10; // VTK_TETRA = 10
      else if(TD == 2)
          return 5; // VTK_TRIANGLE = 5
      else if(TD == 1)
          return 3; // VTK_LINE = 1
      else
          return 0; // VTK_EMPLTY_CELL = 0
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

  void init_top()
  {
      m_face2cell.clear();
      m_cell2face.clear();

      auto NN = number_of_nodes();
      auto NC = number_of_cells();
      m_face2cell.reserve(2*NC); //TODO: 欧拉公式?
      m_cell2face.resize(NC);
      std::map<I, I> idxmap;

      I NF = 0;
      // 偏历所有单元
      for(I i = 0; i < NC; i++)
      {
          for(I j = 0; j < 4; j++)
          {
             auto s = local_face_index(i, j);
             auto it = idxmap.find(s);
             if(it == idxmap.end())
             {
                m_cell2face[i][j] = NF;
                idxmap.insert(std::pair<I, I>(s, NF));
                m_face2cell.push_back(Face2cell{i, i, j, j});
                NF++;
             }
             else
             {
                m_cell2face[i][j] = it->second;
                m_face2cell[it->second][1] = i;
                m_face2cell[it->second][3] = j;
             }
          }
      }
      idxmap.clear();

      m_face.resize(NF);
      for(I i = 0; i < NF; i++)
      {
          auto & c = m_cell[m_face2cell[i][0]];
          auto j = m_face2cell[i][2];
          m_face[i][0] = c[m_localface[j][0]];
          m_face[i][1] = c[m_localface[j][1]];
          m_face[i][2] = c[m_localface[j][2]];
      }

      m_cell2edge.resize(NC);
      I NE = 0;
      for(I i = 0; i < NC; i++)
      {
          for(I j = 0; j < 6; j++)
          { 
              auto & c = m_cell[i];
              auto s = local_edge_index(i, j); 
              auto it = idxmap.find(s);
              if(it == idxmap.end())
              {
                  m_cell2edge[i][j] = NE;
                  idxmap.insert(std::pair<I, I>(s, NE));
                  m_edge.push_back(Edge{c[m_localedge[j][0]], c[m_localedge[j][1]]});
                  NE++; 
              }
             else
             {
                m_cell2edge[i][j] = it->second;
             }
          }
      }
  }

  F cell_quality(const I i)
  {
    auto s = cell_surface_area(i);
    auto d = direction(i, 0);
    auto l = std::sqrt(d.squared_length());
    auto vol = cell_measure(i);
    auto R = l/vol/12.0;
    auto r = 3.0*vol/s;
    return r*3.0/R;
  }

  void cell_quality(std::vector<F> & q)
  {
    auto NC = number_of_cells();
    q.resize(NC);
    for(I i = 0; i < NC; i++)
      q[i] = cell_quality(i);
  }

  F cell_surface_area(const I i)
  {
    auto s = face_measure(m_cell2face[i][0]);
    s += face_measure(m_cell2face[i][1]);
    s += face_measure(m_cell2face[i][2]);
    s += face_measure(m_cell2face[i][3]);
    return s;
  }

  Vector direction(const I i, const I j)
  {
    auto v10 = m_node[m_cell[i][m_index[3*j][0]]] - m_node[m_cell[i][m_index[3*j][1]]];
    auto v20 = m_node[m_cell[i][m_index[3*j][0]]] - m_node[m_cell[i][m_index[3*j][2]]];
    auto v30 = m_node[m_cell[i][m_index[3*j][0]]] - m_node[m_cell[i][m_index[3*j][3]]];

    auto v1 = cross(v20, v30);
    v1 *= v10.squared_length();

    auto v2 = cross(v30, v10);
    v2 *= v20.squared_length();

    auto v3 = cross(v10, v20);
    v3 *= v30.squared_length();

    return v1 + v2 + v3;
  }

  void cell_dihedral_angle(F & max, F & min)
  {
    auto NC = number_of_cells();
    max = 0.0;
    min = 1e+10;
    for(I i=0; i < NC; i++)
    {
      cell_dihedral_angle(i, max, min);
    }
  }

  void cell_dihedral_angle(const I i, F & max, F & min)
  {
    std::array<Vector, 4> ns;
    for(I j = 0; j < 4; j++)
    {
      auto v0 = m_node[m_cell[i][m_localface[j][1]]] - m_node[m_cell[i][m_localface[j][0]]];
      auto v1 = m_node[m_cell[i][m_localface[j][2]]] - m_node[m_cell[i][m_localface[j][0]]];
      ns[j] = cross(v0, v1);
      ns[j] /= std::sqrt(ns[j].squared_length());
    }
    auto PI = GK::pi();
    for(I j = 0; j < 4; j++)
      for(I k = j+1; k < 4; k++)
      {
        F a = PI - std::acos(dot(ns[j], ns[k])); 
        a /= PI;
        a *= 180;
        if( max < a)
        {
            max = a;
        }
        if( min > a)
        {
            min = a;
        }
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
      return m_face.begin();
  }

  FaceIterator face_end()
  {
      return m_face.end();
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

  Node & node(const I i)
  {
      return m_node[i];
  }

  Edge & edge(const I i)
  {
      return m_edge[i];
  }

  Face & face(const I i)
  {
      return m_face[i];
  }

  Cell & cell(const I i)
  {
      return m_cell[i];
  }


  Face2cell & face_to_cell(const I i)
  {
      return m_face2cell[i];
  }

  Cell2face & cell_to_face(const I i)
  {
      return m_cell2face[i];
  }

  Cell2edge & cell_to_edge(const I i)
  {
      return m_cell2edge[i];
  }

  NodeIntData & get_node_int_data()
  {
    return m_nodeintdata;
  }

  NodeDoubleData & get_node_double_data()
  {
    return m_nodedoubledata;
  }
  // 实体测度 

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
      auto & f = m_face[i];
      auto v1 = m_node[f[1]] - m_node[f[0]];
      auto v2 = m_node[f[2]] - m_node[f[0]];
      return 0.5*std::sqrt(cross(v1, v2).squared_length());
  }

  void face_measure(std::vector<F> & measure)
  {
      auto NF = number_of_faces();
      measure.resize(NF);
      for(I i = 0; i < NF; i++)
          measure[i] = face_measure(i);
  }

  F cell_measure(const I i)
  {
      auto & c = m_cell[i];
      auto v01 = m_node[c[1]] - m_node[c[0]];
      auto v02 = m_node[c[2]] - m_node[c[0]];
      auto v03 = m_node[c[3]] - m_node[c[0]];
      return dot(cross(v01, v02), v03)/6.0;
  }

  void cell_measure(std::vector<F> & measure)
  {
      auto NC = number_of_cells();
      measure.resize(NC);
      for(I i = 0; i < NC; i++)
          measure[i] = cell_measure(i);
  }

  // 实体重心
  Node edge_barycenter(const I i)
  {
      auto & e = m_edge[i];
      F x = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
      F y = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
      F z = (m_node[e[0]][2] + m_node[e[1]][2])/2.0;
      return Node(x, y, z);
  }

  void edge_barycenter(const I i, Node & node)
  {
      auto & e = m_edge[i];
      node[0] = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
      node[1] = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
      node[2] = (m_node[e[0]][2] + m_node[e[1]][2])/2.0;
  }

  Node face_barycenter(const I i)
  {
      auto & f = m_face[i];
      F x = (m_node[f[0]][0] + m_node[f[1]][0] + m_node[f[2]][0])/3.0;
      F y = (m_node[f[0]][1] + m_node[f[1]][1] + m_node[f[2]][1])/3.0;
      F z = (m_node[f[0]][2] + m_node[f[1]][2] + m_node[f[2]][2])/3.0;
      return Node(x, y, z);
  }

  void cell_barycenter(const I i, Node & node)
  {
    auto & c = m_cell[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[c[0]][i] + m_node[c[1]][i] + m_node[c[2]][i] + m_node[c[3]][i])/4.0;
  }

  void uniform_refine(const I n=1)
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
          m_cell.resize(8*NC);
          for(I j = 0; j < NC; j++)
          { 
              auto c = m_cell[j]; 
              m_cell[j][0] = c[0];
              m_cell[j][1] = m_cell2edge[j][0] + NN;
              m_cell[j][2] = m_cell2edge[j][1] + NN;
              m_cell[j][3] = m_cell2edge[j][2] + NN;

              m_cell[NC + j][0] = c[1];
              m_cell[NC + j][1] = m_cell2edge[j][3] + NN;
              m_cell[NC + j][2] = m_cell2edge[j][0] + NN;
              m_cell[NC + j][3] = m_cell2edge[j][4] + NN;

              m_cell[2*NC + j][0] = c[2];
              m_cell[2*NC + j][1] = m_cell2edge[j][1] + NN;
              m_cell[2*NC + j][2] = m_cell2edge[j][3] + NN;
              m_cell[2*NC + j][3] = m_cell2edge[j][5] + NN;

              m_cell[3*NC + j][0] = c[3];
              m_cell[3*NC + j][1] = m_cell2edge[j][4] + NN;
              m_cell[3*NC + j][2] = m_cell2edge[j][2] + NN;
              m_cell[3*NC + j][3] = m_cell2edge[j][5] + NN;

              m_cell[3*NC + j][0] = c[3];
              m_cell[3*NC + j][1] = m_cell2edge[j][4] + NN;
              m_cell[3*NC + j][2] = m_cell2edge[j][2] + NN;
              m_cell[3*NC + j][3] = m_cell2edge[j][5] + NN;

              auto v0 = m_node[m_cell2edge[j][0] + NN] - m_node[m_cell2edge[j][5] + NN];
              auto v1 = m_node[m_cell2edge[j][1] + NN] - m_node[m_cell2edge[j][4] + NN]; 
              auto v2 = m_node[m_cell2edge[j][2] + NN] - m_node[m_cell2edge[j][3] + NN]; 

              std::vector<F> v{v0.squared_length(), v1.squared_length(), v2.squared_length()};
              auto idx = std::distance(v.begin(), std::min_element(v.begin(), v.end()));

              m_cell[4*NC + j][0] = m_cell2edge[j][m_refine[idx][0]] + NN; 
              m_cell[4*NC + j][1] = m_cell2edge[j][m_refine[idx][1]] + NN;
              m_cell[4*NC + j][2] = m_cell2edge[j][m_refine[idx][4]] + NN;
              m_cell[4*NC + j][3] = m_cell2edge[j][m_refine[idx][5]] + NN;

              m_cell[5*NC + j][0] = m_cell2edge[j][m_refine[idx][1]] + NN; 
              m_cell[5*NC + j][1] = m_cell2edge[j][m_refine[idx][2]] + NN;
              m_cell[5*NC + j][2] = m_cell2edge[j][m_refine[idx][4]] + NN;
              m_cell[5*NC + j][3] = m_cell2edge[j][m_refine[idx][5]] + NN;

              m_cell[6*NC + j][0] = m_cell2edge[j][m_refine[idx][2]] + NN; 
              m_cell[6*NC + j][1] = m_cell2edge[j][m_refine[idx][3]] + NN;
              m_cell[6*NC + j][2] = m_cell2edge[j][m_refine[idx][4]] + NN;
              m_cell[6*NC + j][3] = m_cell2edge[j][m_refine[idx][5]] + NN;

              m_cell[7*NC + j][0] = m_cell2edge[j][m_refine[idx][3]] + NN; 
              m_cell[7*NC + j][1] = m_cell2edge[j][m_refine[idx][0]] + NN;
              m_cell[7*NC + j][2] = m_cell2edge[j][m_refine[idx][4]] + NN;
              m_cell[7*NC + j][3] = m_cell2edge[j][m_refine[idx][5]] + NN;
          }

          m_edge.clear();
          m_face.clear();
          m_cell2edge.clear();
          m_face2cell.clear();
          m_cell2face.clear();
          init_top(); 
      }
  }

  void is_boundary_face(std::vector<bool> & isBdFace)
  {
      auto NF = number_of_faces();
      isBdFace.resize(NF);
      for(int i = 0; i < NF; i++)
      {
          isBdFace[i] = false;
          if(face_to_cell(i)[0]==face_to_cell(i)[1])
          {
              isBdFace[i] = true;
          }
      }
  }

  void is_boundary_node(std::vector<bool> & isBdNode)
  {
      auto NN = number_of_nodes();
      auto NF = number_of_faces();
      isBdNode.resize(NN);
      for(int i = 0; i < NF; i++)
      {
          if(face_to_cell(i)[0]==face_to_cell(i)[1])
          {
              isBdNode[face(i)[0]] = true;
              isBdNode[face(i)[1]] = true;
              isBdNode[face(i)[2]] = true;
          }
      }
  }

  void cell_to_node(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NN = number_of_nodes();
      auto nn = number_of_vertices_of_each_cell();

      top.init(3, 0, NC, NN);

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
          nei[nn*i + 3] = m_cell[i][3];
      }
  }

  void cell_to_cell(Toplogy & top)
  {

      auto NC = number_of_cells();
      auto NF = number_of_faces();
      top.init(3, 3, NC, NC);
      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      for(I i=0; i < NF; i++)
      {
          loc[m_face2cell[i][0]+1] += 1;
          if(m_face2cell[i][0] != m_face2cell[i][1])
          {
              loc[m_face2cell[i][1]+1] += 1;
          }
      }
      for(I i=0; i < NC; i++)
      {
          loc[i+1] += loc[i];
      }

      auto & nei = top.neighbors();
      nei.resize(loc[NC]);
      std::vector<I> start = loc;
      for(I i = 0; i < NF; i++)
      {
          nei[start[m_face2cell[i][0]]++] = m_face2cell[i][1];
          if(m_face2cell[i][0] != m_face2cell[i][1])
          {
              nei[start[m_face2cell[i][1]]++] = m_face2cell[i][0];
          }
      }
  }

  void node_to_node(Toplogy & top)
  {
      auto NN = number_of_nodes();
      auto NE = number_of_edges();
      top.init(3, 3, NN, NN);
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
    auto & locid = top.local_indices();
    nei.resize(loc[NN]);
    locid.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NC; i++)
    {
      for(auto j = 0; j < 4; j++)
      {
        auto v = m_cell[i][j];
        locid[start[v]] = j;
        nei[start[v]++] = i;
      }
    }
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      print_entities(m_node);

      std::cout << "Edges:" << std::endl;
      print_entities(m_edge);

      std::cout << "Faces:" << std::endl;
      print_entities(m_face);

      std::cout << "Cells:" << std::endl;
      print_entities(m_cell);

      std::cout << "Face2cell:" << std::endl;
      print_entities(m_face2cell);

      std::cout << "Cell2face:" << std::endl;
      print_entities(m_cell2face);
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
     *  计算第 i 个 cell 的第 j 个 face 的全局唯一整数编号
     */
    I local_face_index(I i, I j)
    {
        I f[3] = {
            m_cell[i][m_localface[j][0]], 
            m_cell[i][m_localface[j][1]], 
            m_cell[i][m_localface[j][2]]
        };
        std::sort(f, f+3);
        return  f[0] + f[1]*(f[1]+1)/2 + f[2]*(f[2]+1)*(f[2]+2)/6;
    }

    /*
     *
     * Notes
     * -----
     *  计算第 i 个 cell 的 j 条 edge 全局唯一整数编号
     */
    I local_edge_index(I i, I j)
    {
        I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }

private:
    std::vector<Node> m_node;
    std::vector<Cell> m_cell; 
    std::vector<Edge> m_edge;
    std::vector<Face> m_face;
    std::vector<Face2cell> m_face2cell;
    std::vector<Cell2face> m_cell2face;
    std::vector<Cell2edge> m_cell2edge;
    NodeIntData m_nodeintdata;
    NodeDoubleData m_nodedoubledata;
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_localedge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_localface[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_localface2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_refine[3][6] = {
    {1, 3, 4, 2, 5, 0}, {0, 2, 5, 3, 4, 1}, {0, 4, 5, 1, 3, 2}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_index[12][4] = {
    {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
    {1, 2, 0, 3}, {1, 0, 3, 2}, {1, 3, 2, 0},
    {2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1},
    {3, 0, 2, 1}, {3, 2, 1, 0}, {3, 1, 0, 2}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_vtkidx[4] = {0, 1, 2, 3};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::m_num[4][4] = {
    {0, 1, 2, 3}, {1, 0, 3, 2}, {2, 0, 1, 3}, {3, 0, 2, 1}
};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of TetrahedronMesh_h

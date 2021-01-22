#ifndef TetrahedronMesh_h
#define TetrahedronMesh_h

#include <vector>
#include <array>
#include <map>

#include "MeshToplogy.h"

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
template<typename GK, typename Node, typename Vector>
class TetrahedronMesh 
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;

    // 拓扑
    typedef MeshToplogy<I> Toplogy;

    // 实体类型
    typedef typename std::array<I, 4> Cell;
    typedef typename std::array<I, 2> Edge;
    typedef typename std::array<I, 3> Face;

    // 实体关系类型
    typedef typename std::array<I, 4> Face2cell;
    typedef typename std::array<I, 4> Cell2Face;
    typedef typename std::array<I, 6> Cell2edge;

    // 迭代子类型
    typedef typename std::vector<Node>::iterator Node_iterator;
    typedef typename std::vector<Edge>::iterator Edge_iterator;
    typedef typename std::vector<Face>::iterator Face_iterator;
    typedef typename std::vector<Cell>::iterator Cell_iterator;

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
        return m_edge.size();
    }

    static int geo_dimension()
    {
        return 3;
    }

    static int top_dimension()
    {
        return 3;
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


    Node_iterator node_begin()
    {
        return m_node.begin();
    }

    Node_iterator node_end()
    {
        return m_node.end();
    }

    Edge_iterator edge_begin()
    {
        return m_edge.begin();
    }

    Edge_iterator edge_end()
    {
        return m_edge.end();
    }

    Face_iterator face_begin()
    {
        return m_face.begin();
    }

    Face_iterator face_end()
    {
        return m_face.end();
    }

    Cell_iterator cell_begin()
    {
        return m_cell.begin();
    }

    Cell_iterator cell_end()
    {
        return m_cell.end();
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
        auto v1 = m_node[c[1]] - m_node[c[0]];
        auto v2 = m_node[c[2]] - m_node[c[0]];
        auto v3 = m_node[c[3]] - m_node[c[0]];
        return dot(cross(v1, v2), v3)/6.0;
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
        return Node(x, y);
    }

    void edge_barycenter(const I i, Node & node)
    {
        auto & e = m_edge[i];
        node[0] = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
        node[1] = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
    }

    Node face_barycenter(const I i)
    {
        auto & f = m_face[i];
        F x = (m_node[f[0]][0] + m_node[f[1]][0] + )/2.0;
        F y = (m_node[f[0]][1] + m_node[e[1]][1])/2.0;
        F z = 
        return Node(x, y);
    }

    void edge_barycenter(const I i, Node & node)
    {
        auto & e = m_edge[i];
        node[0] = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
        node[1] = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
    }

    void uniform_refine(const int n=1)
    {
        for(I i=0; i < n; i++)
        {
            m_edge.clear();
            m_face.clear();
            m_cell2edge.clear();
            m_face2cell.clear();
            init_top(); 
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
     *  计算第 i 个 cell 的第 j-th face 的全局唯一整数编号
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
     *  计算第 i 个 cell 的 j-th edge 全局唯一整数编号
     */
    I local_edge_index(I i, I j)
    {
        I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }

private:
    static int m_localedge[6][2];
    static int m_localface[4][3];
    static int m_localface2edge[4][3];
    std::vector<Node> m_node;
    std::vector<Cell> m_cell; 
    std::vector<Edge> m_edge;
    std::vector<Face> m_face;
    std::vector<Edge2cell> m_face2cell;
    std::vector<Cell2face> m_cell2face;
    std::vector<Cell2edge> m_cell2edge;
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::localedge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::localface[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

template<typename GK, typename Node, typename Vector>
int TetrahedronMesh<GK, Node, Vector>::localface2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of TetrahedronMesh_h
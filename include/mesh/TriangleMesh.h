#ifndef TriangleMesh_h
#define TriangleMesh_h

#include <vector>
#include <array>
#include <map>

#include "MeshToplogy.h"

namespace WHYSC {
namespace Mesh {

/*
 * 
 *
 *
 */
template<typename GK, typename Node, typename Vector>
class TriangleMesh 
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;

    typedef MeshToplogy<I> Toplogy;
    typedef typename std::array<I, 3> Cell;
    typedef typename std::array<I, 2> Edge;
    typedef Edge Face;
    typedef typename std::array<I, 4> Edge2cell;
    typedef typename std::array<I, 3> Cell2edge;

public:
    TriangleMesh()
    {
    }

    void insert(Node & node)
    {
        m_node.push_back(node);
    }

    void insert(Cell & cell)
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

    int geo_dimension()
    {
        return Node::dimension();
    }

    static int top_dimension()
    {
        return 2;
    }

    void construct_top()
    {
        I NE = 0;
        std::map<I, I> idxmap;
        m_cell2edge.resize(m_cell.size());

        // 偏历所有单元
        for(I i = 0; i < m_cell.size(); i++)
        {
            for(I j=0; j < 3; j++)
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
    }

    /*
     *
     * Notes
     *  计算第 i 个单元的第 j 条边全局唯一的一个整数索引
     */
    I local_edge_index(I i, I j)
    {
        I e[2] = {m_cell[i][edge[j][0]], m_cell[i][edge[j][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }
private:
    static int ccw[3];
    static int cw[3];
    static int edge[3][2];
    static int face[3][2];

    std::vector<Node> m_node;
    std::vector<Cell> m_cell; 
    std::vector<Edge> m_edge;
    std::vector<Edge2cell> m_edge2cell;
    std::vector<Cell2edge> m_cell2edge;
};

template<typename GK, typename Node, typename Vector>
int TriangleMesh<GK, Node, Vector>::cw[3] = {2, 0, 1};

template<typename GK, typename Node, typename Vector>
int TriangleMesh<GK, Node, Vector>::ccw[3] = {1, 2, 0};


template<typename GK, typename Node, typename Vector>
int TriangleMesh<GK, Node, Vector>::edge[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};

template<typename GK, typename Node, typename Vector>
int TriangleMesh<GK, Node, Vector>::face[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of TriangleMesh_h

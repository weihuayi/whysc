#ifndef Hexahedron_h
#define Hexahedron_h

#include <algorithm>
#include <utility>
#include <initializer_list>
#include "Cell_type.h"

namespace WHYSC {
namespace Mesh {

template<typename I>
struct Hexahedron
{
    struct Face 
    {
        Face(const std::initializer_list<I> &l)
        {
            std::copy_n(l.begin(), 2, m_cell);
            std::copy_n(l.begin()+2, 2, m_index);
        }

        I & cell(const int i)
        {
            return m_cell[i];
        }

        const I & cell(const int i) const
        {
            return m_cell[i];
        }

        I & index(const int i)
        {
            return m_index[i];
        }

        const I & index(const int i) const
        {
            return m_index[i];
        }
    private:
        I m_cell[2];
        I m_index[2];
    };

    struct Edge
    {
        Edge(const std::initializer_list<I> &l)
        {
            std::copy_n(l.begin(), 1, m_cell);
            std::copy_n(l.begin()+1, 1, m_index);
        }

        I & localnode(int i)
        {
            return Hexahedron::localedge[m_index[0]][i];
        }

        I & cell(const int i)
        {
            return m_cell[i];
        }

        const I & cell(const int i) const
        {
            return m_cell[i];
        }

        I & index(const int i)
        {
            return m_index[i];
        }

        const I & index(const int i) const
        {
            return m_index[i];
        }

    private:
        I m_cell[1];
        I m_index[1];
    };

    static int NV[4];
    static int ND[4];
    static Cell_type type;

    static int localedge[12][2];
    static int localface[6][4];
    static int localface2edge[6][4];

    Hexahedron(const std::initializer_list<I> &l)
    { 
        std::copy_n(l.begin(), 8, m_node);
    }

    I & node(const int i)
    {
        return m_node[i];
    }

    const I & node(const int i) const
    {
        return m_node[i];
    }

    I & edge(const int i)
    {
        return m_edge[i];
    }

    const I & edge(const int i) const
    {
        return m_edge[i];
    }

    I & face(const int i)
    {
        return m_face[i];
    }

    const I & face(const int i) const
    {
        return m_face[i];
    }


    static int dimension() {return 3;}

    I local_face_index(int i)
    {
        I f[4] = {m_node[localface[i][0]], m_node[localface[i][1]], m_node[localface[i][2]],
            m_node[localface[i][3]]};
        std::sort(f, f+4);
        return f[0] + f[1]*(f[1]+1)/2 + f[2]*(f[2] + 1)*(f[2] + 2)/6 + f[3]*(f[3]+1)*(f[3]+2)*(f[3]+3)/24;
    }

    I local_edge_index(int i)
    {
        I e[2] = {m_node[localedge[i][0]], m_node[localedge[i][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }

private:
    I m_node[8];
    I m_edge[12];
    I m_face[6];
};

template<typename I>
int Hexahedron<I>::NV[4] = {1, 2, 4, 8};

template<typename I>
int Hexahedron<I>::ND[4] = {8, 12, 6, 1};

template<typename I>
Cell_type Hexahedron<I>::type = 12; // VTK_TETRA

template<typename I>
int Hexahedron<I>::localedge[12][2] = {
    {0, 1}, {0, 2}, {0, 4}, {1, 3}, 
    {1, 5}, {2, 3}, {2, 6}, {3, 7},
    {4, 5}, {4, 6}, {5, 7}, {6, 7}
};

template<typename I>
int Hexahedron<I>::localface[6][4] = {
    {0, 2, 6, 4}, {1, 5, 7, 3},
    {0, 1, 3, 2}, {4, 6, 7, 5},  
    {0, 4, 5, 1}, {6, 2, 3, 7}
};

template<typename I>
int Hexahedron<I>::localface2edge[6][4] = {
    {1,  6, 9, 2}, {4, 10, 7, 3},
    {0,  3, 5, 1}, {9, 11, 10, 8},
    {2,  8, 4, 0}, {6, 5, 7,  11}
};
} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of cells_h

#ifndef cells_h
#define cells_h

#include <algorithm>
#include <utility>
#include <initializer_list>
#include "Cell_type.h"

namespace WHYSC {
namespace Mesh {

template<typename I>
struct Tetrahedron
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
            std::copy_n(l.begin(), 2, m_node);
        }

        I & cell(const int i)
        {
            return m_cell;
        }

        const I & cell(const int i) const
        {
            return m_cell;
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

    static int localedge[6][2];
    static int localface[4][3];
    static int localface2edge[4][3];

    Tetrahedron(const std::initializer_list<I> &l)
    { 
        std::copy_n(l.begin(), 4, m_node);
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
        I f[3] = {m_node[localface[i][0]], m_node[localface[i][1]], m_node[localface[i][1]]};
        std::sort(f, f+3, std::greater<I>());
        return f[0] + f[1]*(f[1]+1)/2 + f[2]*(f[2] + 1)*(f[2] + 2)/6;
    }

    I local_edge(int i)
    {
        I e[2] = {m_node[localface[i][0]], m_node[localface[i][1]]};
        std::sort(e, e+2, std::greater<I>());
        return  e[0] + e[1]*(e[1]+1)/2;
    }

private:
    I m_node[4];
    I m_edge[6];
    I m_face[4];
};

template<typename I>
int Tetrahedron<I>::NV[4] = {1, 2, 3, 4};

template<typename I>
int Tetrahedron<I>::ND[4] = {4, 6, 4, 1};

template<typename I>
Cell_type Tetrahedron<I>::type = 10; // VTK_TETRA

template<typename I>
int Tetrahedron<I>::localedge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

template<typename I>
int Tetrahedron<I>::localface[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

template<typename I>
int Tetrahedron<I>::localface2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of cells_h

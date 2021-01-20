#ifndef Triangle_h
#define Triangle_h

#include <initializer_list>
#include "Cell_type.h"

namespace WHYSC {

namespace Mesh {

template<typename Point>
class Node: public Point
{

pivate:
    int m_id;
};


template<typename Node>
class Triangle
{
    struct Face 
    {
        Face(const std::initializer_list<int> &l)
        { 
            std::copy_n(l.begin(), 2, m_cell);
            std::copy_n(l.begin()+2, 2, m_index);
        }

        I & node(const int i, Triangle<I> * cells)
        {
            return cells[m_cell[0]].node(localface[m_index[0]][i]);
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

    typedef Face Edge;


    static int NV[3]; 
    static int ND[3];
    static Cell_type type;
    static int localedge[3][2];
    static int localface[3][2];
    static int localface2edge[3][1];

    Triangle(const std::initializer_list<I> &l)
    { 
        std::copy_n(l.begin(), 3, m_node);
    }

    I & node(const int i)
    {
        return m_node[i];
    }

    const I & node(const int i) const
    {
        return m_node[i];
    }

    I & face(const int i)
    {
        return m_face[i];
    }

    const I & face(const int i) const
    {
        return m_face[i];
    }

    I & edge(const int i)
    {
        return m_face[i];
    }

    const I & edge(const int i) const
    {
        return m_face[i];
    }

    static int dimension() {return 2;}

    I local_face_index(int i)
    {
        I f[2] = {m_node[localface[i][0]], m_node[localface[i][1]]};
        std::sort(f, f+2);
        return  f[0] + f[1]*(f[1]+1)/2;
    }

    I local_edge_index(int i)
    {
        I e[2] = {m_node[localface[i][0]], m_node[localface[i][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }

private:
    I m_node[3];
    I m_face[3];
};


template<typename I>
int Triangle<I>::NV[3] = {1, 2, 3};

template<typename I>
int Triangle<I>::ND[3] = {3, 3, 1};

template<typename I>
Cell_type Triangle<I>::type = 5; // VTK_TRIANGLE


template<typename I>
int Triangle<I>::localedge[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};

template<typename I>
int Triangle<I>::localface[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};

template<typename I>
int Triangle<I>::localface2edge[3][1] = {
    {0}, {1}, {2}
};

template<typename I>
std::ostream& operator << (std::ostream & os, const Triangle<I> & cell)
{

    auto TD = Triangle<I>::dimension();
    for(auto i = 0; i < Triangle<I>::NV[TD]; i++)
    {
        os << cell.node(i) << " ";
    }
    return os;
}

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Triangle_h

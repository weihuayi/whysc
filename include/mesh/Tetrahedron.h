#ifndef cells_h
#define cells_h

#include <array>
#include <algorithm>
#include <utility>
#include "Cell_type.h"

namespace WHYSC {
namespace Mesh {

template<typename I>
struct Tetrahedron: public std::array<I, 4>
{
    typedef typename std::array<I, 4> Base;
    using Base::data;

    typedef std::array<I, 2> Edge;
    struct Face: public std::array<I, 3> 
    {
        I face2cell[4];
    };

    static int dim; // the dimension of cell
    static int NV[4];
    static int ND[4];
    static Cell_type type;

    static int localedge[6][2];
    static int localface[4][3];
    static int localface2edge[4][3];


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

    static int dimension() {return 3;}

    Tetrahedron(const std::initializer_list<I> &l)
    { 
        std::copy_n(l.begin(), 2, data());
    }

    I get_local_face(int i, Face & f)
    {
        f[0] = *this[localface[i][0]];
        f[1] = *this[localface[i][1]];
        f[2] = *this[localface[i][2]];

        Face f0 = f;
        std::sort(f0.begin(), f0.end(), std::greater<I>());
        return f0[0] + f0[1]*(f0[1]+1)/2 + f0[2]*(f0[2] + 1)*(f0[2] + 2)/6;
    }

    I get_local_edge(int i, Edge & e)
    {
        e[0] = *this[localedge[i][0]];
        e[1] = *this[localedge[i][1]];
        I s;
        if (e[0] < e[1])
            s = e[0] + e[1]*(e[1]+1)/2;
        else
            s = e[1] + e[0]*(e[0]+1)/2;
        return s;
    }

private:
    I m_edge[6];
    I m_face[4];
};

template<typename I>
int Tetrahedron<I>::NV[4] = {1, 2, 3, 4};

template<typename I>
int Tetrahedron<I>::ND[4] = {4, 6, 4, 1};

template<typename I>
Cell_type Tetrahedron<I>::type = TETRA;

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

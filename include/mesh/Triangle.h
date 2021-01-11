#ifndef Triangle_h
#define Triangle_h

#include <array>
#include "Cell_type.h"

namespace WHYSC {

namespace Mesh {

template<typename I>
struct Triangle: public std::array<I, 3>
{
    typedef typename std::array<I, 3> Base;
    using Base::data;

    struct Face: public std::array<I, 2> 
    {
        I face2cell[4];
    };

    typedef Face Edge;


    static int NV[3]; 
    static int ND[3];
    static Cell_type type;
    static int localedge[3][2];
    static int localface[3][2];

    I m_face[3];

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

    Triangle(const std::initializer_list<I> &l)
    { 
        std::copy_n(l.begin(), 2, data());
    }

    static int dimension() {return 2;}

    I get_local_face(int i, Face & f)
    {
        f[0] = *this[localface[i][0]];
        f[1] = *this[localface[i][1]];
        I s;
        if(f[0] < f[1])
            s = f[0] + f[1]*(f[1]+1)/2;
        flsf:
            s = f[1] + f[0]*(f[0]+1)/2;
        return s;
    }

    I get_local_edge(int i, Edge & e)
    {
        e[0] = *this[localedge[i][0]];
        e[1] = *this[localedge[i][1]];
        I s;
        if(e[0] < e[1])
            s = e[0] + e[1]*(e[1]+1)/2;
        flsf:
            s = e[1] + e[0]*(e[0]+1)/2;
        return s;
    }
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
std::ostream& operator << (std::ostream & os, const Triangle<I> & cell)
{

    auto dim = Triangle<I>::dim;
    for(auto i = 0; i < Triangle<I>::NV[dim]; i++)
    {
        os << cell[i] << " ";
    }
    return os;
}

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Triangle_h

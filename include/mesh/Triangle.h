#ifndef Triangle_h
#define Triangle_h

#include <array>
#include "Cell_type.h"

namespace WHYSC {

namespace Mesh {

template<typename I>
struct Triangle: public std::array<I, 3>
{
    typedef std::array<I, 2> Edge;
    typedef Edge Face;

    static int dim; // the toplogy dimension of cell
    static int NV[3]; 
    static int ND[3];
    static CellType type;
    static int edge[3][2];
};

template<typename I>
int Triangle<I>::dim = 2;


template<typename I>
int Triangle<I>::NV[3] = {1, 2, 3};

template<typename I>
int Triangle<I>::ND[3] = {3, 3, 1};

template<typename I>
CellType Triangle<I>::type = 5; // VTK_TRIANGLE


template<typename I>
int Triangle::edge[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};


template<typename I>
std::ostream& operator << (std::ostream & os, const Triangle<I> & cell)
{

    auto dim = Triangle::dim;
    for(auto i = 0; i < Triangle::NV[dim]; i++)
    {
        os << cell[i] << " ";
    }
    return os;
}

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Triangle_h

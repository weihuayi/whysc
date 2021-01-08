#ifndef cells_h
#define cells_h

#include <array>
#include <utility>
#include "Cell_type.h"

namespace WHYSC {
namespace Mesh {

template<typename I>
struct Tetrahedron
{
    typedef std::array<I, 2> Edge;
    typedef std::array<I, 3> Face;
    typedef std::pair<int, Face> IFace;
    typedef std::array<int, 4> Cell2Cell;
    typedef std::array<int, 4> Cell2Face;
    typedef std::array<int, 6> Cell2Edge;
    typedef std::array<int, 4> Face2Cell;

    static int dim; // the dimension of cell
    static int NV[4];
    static int ND[4];
    static Cell_type type;

    static int edge[6][2];
    static int face[4][3];
    static int face2edge[4][3];
    
    int _data[4];


    int & operator[](const int i) 
    {
        return _data[i];
    }

    const int & operator[](const int i) const
    {
        return _data[i];
    }
};

template<typename I>
int Tetrahedron<I>::dim = 3;

template<typename I>
int Tetrahedron<I>::NV[4] = {1, 2, 3, 4};

template<typename I>
int Tetrahedron<I>::ND[4] = {4, 6, 4, 1};

template<typename I>
Cell_type Tetrahedron<I>::type = TETRA;

template<typename I>
int Tetrahedron<I>::edge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

template<typename I>
int Tetrahedron<I>::face[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

template<typename I>
int Tetrahedron<I>::face2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

//std::ostream& operator << (std::ostream & os, const Tetrahedron & cell)
//{
//
//    for(auto i = 0; i < Cell::V; i++)
//    {
//        os << cell[i] << " ";
//    }
//    return os;
//}

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of cells_h

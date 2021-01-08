#ifndef cells_h
#define cells_h

#include <array>
#include <utility>
#include "CellType.h"

namespace WHYSC {
namespace Mesh {

struct Tetrahedron
{
    typedef std::array<int, 2> Edge;
    typedef std::array<int, 3> Face;
    typedef std::pair<int, Face> IFace;
    typedef std::array<int, 4> Cell2cell;
    typedef std::array<int, 4> Cell2face;
    typedef std::array<int, 6> Cell2edge;
    typedef std::array<int, 4> Face2cell;

    static int dim; // the dimension of cell
    static int NV[4];
    static int ND[4];
    static CellType type;

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

int Tetrahedron::dim = 3;
int Tetrahedron::NV[4] = {1, 2, 3, 4};
int Tetrahedron::ND[4] = {4, 6, 4, 1};
CellType Tetrahedron::type = TETRA;

int Tetrahedron::edge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

int Tetrahedron::face[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

int Tetrahedron::face2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

std::ostream& operator << (std::ostream & os, const Tetrahedron & cell)
{

    for(auto i = 0; i < Cell::V; i++)
    {
        os << cell[i] << " ";
    }
    return os;
}

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of cells_h

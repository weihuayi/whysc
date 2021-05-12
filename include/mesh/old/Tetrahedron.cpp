#include "Tetrahedron.h"

namespace  iMath {
namespace GeometryObject {
int Tetrahedron::face[4][3]={{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
int Tetrahedron::edge[6][2]= {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
int Tetrahedron::edge_idx[4][4]={
       {0, 0, 1, 2},
       {0, 0, 3, 4},
       {1, 3, 0, 5},
       {2, 4, 5, 0}};
int Tetrahedron::index[12][4]= {
       {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
       {1, 2, 0, 3}, {1, 0, 3, 2}, {1, 3, 2, 0},
       {2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1},
       {3, 0, 2, 1}, {3, 2, 1, 0}, {3, 1, 0, 2}};
}
}

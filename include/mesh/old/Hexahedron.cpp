#include "Hexahedron.h"

namespace iMath {
namespace GeometryObject {

int Hexahedron::face[6][4]={
        0, 3, 2, 1, 4, 5, 6, 7,
        0, 4, 7, 3, 1, 2, 6, 5,
        0, 1, 5, 4, 3, 7, 6, 2};

int Hexahedron::corner[8][3] = {
        1, 3, 4, 2, 0, 5, 3, 1, 6, 0, 2, 7,
        7, 5, 0, 4, 6, 1, 5, 7, 2, 6, 4, 3};
}

}

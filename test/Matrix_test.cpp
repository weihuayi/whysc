
#include <iostream>

#include "algebra/Matrix.h"
#include "algebra/Vector.h"

typedef WHYSC::AlgebraObject::Matrix<double, int> Matrix;
typedef WHYSC::AlgebraObject::Vector<double, int> Vector;

int main(int argc, char **argv)
{

    Matrix M = {{1, 2, 4}, {3, 4, 5}, {6, 7, 8}};
    std::cout << M << std::endl;
    std::cout << "M(1, 2)=" << M(1, 2) << std::endl;
    std::cout << "M[1][2]=" << M[1][2] << std::endl;
    Vector v = {1, 2, 3};
    std::cout << v.norm() << std::endl;


    Matrix M0(10, 10);
    M0.fill_diag(2, 0); // 主对角线填充为 2
    M0.fill_diag(-1, -1);
    M0.fill_diag(-1, 1);

    std::cout << M0 << std::endl;

    return 0;
}

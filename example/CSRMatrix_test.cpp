#include <iostream>

#include "algebra/Algebra_kernel.h"

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;


int main(int argc, char **argv)
{
    Matrix M{{4, 0, 6}, {0, 0, 0}, {0, 0, 5}};
    std::cout << "M:\n" << M << std::endl;

    Vector v0{1, 1, 1};

    CSRMatrix SM(M);
    std::cout << SM << std::endl;
    std::cout << SM.format << std::endl;
    std::cout << SM*v0 << std::endl;

    double data[3] = {4, 6, 5};
    int row[3] = {0, 0, 2};
    int col[3] = {0, 2, 2};
    CSRMatrix SM1(3, 3, 3, row, col, data);

    std::cout << SM1 << std::endl;

    return 0;
}

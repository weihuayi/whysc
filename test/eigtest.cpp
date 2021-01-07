#include <iostream>

#include "algebra/Algebra_kernel.h"

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;
typedef typename Kernel::Eigenvalue Eigenvalue;


int main(int argc, char **argv)
{
    Matrix M{{4, 0, 6}, {0, 0, 0}, {0, 0, 5}};
    auto G = Eigenvalue(M);
    float d = *G.eig;
    std::cout << d << std::endl;
} 

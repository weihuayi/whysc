#include <iostream>

#include "algebra/Algebra_kernel.h"

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;

int main(int argc, char **argv)
{
    Matrix M = {{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}};

    Matrix Q(3, 3);
    Matrix R(3, 3);

    Kernel::qr_gs(M, Q, R);

    std::cout << "M:\n" << M << std::endl;
    std::cout << "Q:\n" << Q << std::endl;
    std::cout << "R:\n" << R << std::endl;
    std::cout << "QR:\n" << Q*R << std::endl;
    return 0;
}

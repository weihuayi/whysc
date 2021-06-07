
#include <iostream>

#include "algebra/Algebra_kernel.h"

#include "algebra/linalg.h"

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;

int main(int argc, char **argv)
{
    Matrix M = {{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}};

    Vector x = {12, -51, 4};
    Vector v(3, 0.0);

    double beta ;
    Kernel::householder(x, v, beta);

    std::cout << "beta: " << beta << std::endl;

    std::cout << "householder vector: " << v << std::endl;

    return 0;
}

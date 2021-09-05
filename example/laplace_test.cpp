#include <iostream>
#include <cmath>

#include "algebra/Algebra_kernel.h"
#include "Constants.h"

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;

int main(int argc, char **argv)
{
    int n = 40;
    double h = 1.0/(n+1);
    CSRMatrix A;
    Kernel::laplace_1d(n, A);

    Vector b(n);
    Vector x(n);
    Vector u(n);

    double pi = WHYSC::Constants::pi;
    for(int i = 0; i < n; i++)
    {
        double x = (i+1)*h;
        u[i] = std::sin(pi*x);
        b[i] = h*h*pi*pi*u[i];
    }

    std::cout << b << std::endl;
    A.jacobi(b, x, 10000);

    //Vector r = b - A*x;
    //std::cout << r.norm() << std::endl;

    double e = (u - x).maxnorm();
    std::cout << "error: " << e << std::endl;

    return 0;
}

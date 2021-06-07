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
    int n = 64;
    int N = n*n;
    double h = 1.0/(n+1);
    CSRMatrix A;
    Kernel::laplace_2d(N, A);

    Vector b(N);
    Vector x(N);
    Vector u(N);

    double pi = WHYSC::Constants::pi;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
        {
            double x = (i+1)*h;
            double y = (j+1)*h;
            u[i] = std::sin(pi*x)*std::sin(pi*y);
            b[i] = 2*h*h*pi*pi*u[i];
        }

    std::cout << b << std::endl;
    A.jacobi(b, x, 10000);

    double e = (u - x).maxnorm();
    std::cout << "error: " << e << std::endl;

    return 0;
}

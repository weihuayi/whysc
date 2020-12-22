#include <iostream>
#include <cstdlib>
#include <ctime>

#include "algebra/Algebra_kernel.h"

// rand(), 0 to RAND_MAX 
// srand(), 设定随机数种子

typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename Kernel::Matrix Matrix;

int main()
{
    srand((unsigned) time(NULL));
    std::cout << "RAND_MAX: " <<  RAND_MAX << std::endl;
    for(int i = 0; i < 5; i++)
    {
        std::cout << rand() << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < 5; i++)
    {
        std::cout << (double) rand()/RAND_MAX << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < 5; i++)
    {
        std::cout << rand()%10 + 1 << std::endl;
    }

    Matrix M(10, 10);
    Kernel::randn_matrix(M, -10, 10);
    std::cout << "M: " << M << std::endl;
}


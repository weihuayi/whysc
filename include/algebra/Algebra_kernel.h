#ifndef AlgebraKernel_h
#define AlgebraKernel_h

#include "Matrix.h"
#include "CSRMatrix.h"
#include "Vector.h"
#include "linalg.h"
#include "MatrixFactory.h"

namespace WHYSC {

template<typename F=double, typename I=int>
class Algebra_kernel
{
public:
    typedef typename AlgebraObject::Matrix<F, I> Matrix;
    typedef typename AlgebraObject::Vector<F, I> Vector;

    typedef typename AlgebraObject::CSRMatrix<F, I> CSRMatrix;
    //typedef AlgebraAlgrithom::lu<Matrix> lu;
public:
    static void lu(Matrix & A, Matrix & L, Matrix & U)
    {
        return AlgebraAlgrithom::lu<Matrix>(A, L, U);
    }

    static void qr_gs(Matrix & A, Matrix & Q, Matrix & R)
    {
        return AlgebraAlgrithom::qr_gs<Matrix>(A, Q, R);
    }

    static void householder(Vector & x, Vector & v, F & beta)
    {
        return AlgebraAlgrithom::householder(x, v, beta);
    }

    static void givens(F x0, F x1, F & c, F & s)
    {
        return AlgebraAlgrithom::givens(x0, x1, c, s);
    }

    static void randn_matrix(Matrix & M, 
            int low=0, 
            int upper=10,
            int seed=0)
    {
        return AlgebraObject::randn_matrix(M, low, upper, seed);
    }

    static void laplace_1d(const I n, CSRMatrix & m)
    {
        return AlgebraObject::laplace_1d(n, m);
    }

    static void laplace_2d(const I N, CSRMatrix & m)
    {
        return AlgebraObject::laplace_2d(N, m);
    }

    static void laplace_3d(const I M, CSRMatrix & m)
    {
        return AlgebraObject::laplace_3d(M, m);
    }
};

} // end of namespace WHYSC
#endif // end of AlgebraKernel_h

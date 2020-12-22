#ifndef linalg_h
#define linalg_h


#include <cstdlib>
#include <ctime>
#include <iostream>
#include <initializer_list>

namespace WHYSC {
namespace AlgebraAlgrithom {

template <typename T> int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template<typename Vector, typename F>
void householder(Vector & x, Vector & v, F & beta)
{
    typedef typename Vector::Int  I;
    F sigma = 0.0;

    I n = x.size;
    v[0] = x[0];
    for(I i = 1; i < n; i++)
    {
        v[i] = x[i];
        sigma += x[i]*x[i];
    }

    if(sigma == 0.0)
    {
        if(x[0] < 0)
        {
            v[0] = 2*x[0];
            beta = 2/(v[0]*v[0]);
        }
        else
        {
            v[0] = 0.0;
            beta = 0.0;
        }

    }
    else
    {
        F alpha = std::sqrt(x[0]*x[0] + sigma);
        if(x[0] < 0.0)
        {
            v[0] -= alpha;
        }
        else
        {
            v[0] = -sigma/(x[0] + alpha);
        }
        beta = 2.0/(v[0]*v[0] + sigma);
    }

    return;
}

template<typename F>
void givens(F x0, F x1, F & c, F & s)
{
    if(x1 == 0.0)
    {
        if( x0 >= 0.0 )
        {
            c = 1.0;
        }
        else
        {
            c = -1.0;
        }
        s = 0.0;
    }
    else
    {
        if(std::abs(x1) > std::abs(x0))
        {
            F tau = x0/x1;
            s = sign(x1)/std::sqrt(1 + tau*tau);
            c = s*tau;
        }
        else
        {
            F tau = x1/x0;
            c = sign(x0)/std::sqrt( 1 + tau*tau);
            s = c*tau;
        }
    }
}


template<typename Matrix>
void lu(Matrix & A, Matrix & L, Matrix & U)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			L[i][k] = A[i][k]/A[k][k];
		}

		for(I j=k; j<n; j++)
		{
			U[k][j] = A[k][j];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= L[i][k]*U[k][j];
			}
		}

	}

    for(I i=0;i<n;i++)
        L[i][i] =1 ;

    return;
}

template<typename Matrix>
void lu_kij(Matrix & A)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			A[i][k] = A[i][k]/A[k][k];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= A[i][k]*A[k][j];
			}
		}

	}

    return;
}

template<typename Matrix>
void lu_ikj(Matrix & A)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I i = 0;i < n; i++)
    	{
   		for(I k = 0;k < i; k++)
   		{
   			A[i][k]=A[i][k]/A[k][k];
   			for(I j = k + 1; j < n; j++)
   			{
   				A[i][j] -= A[i][k]*A[k][j];
   			}
   		}
    	}
    return;
}

template<typename Matrix>
void lu_kij(Matrix & A, Matrix & L, Matrix & U)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			L[i][k] = A[i][k]/A[k][k];
		}

		for(I j=k; j<n; j++)
		{
			U[k][j] = A[k][j];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= L[i][k]*U[k][j];
			}
		}

	}

    for(I i=0;i<n;i++)
        L[i][i] =1 ;

    return;
}


template<typename Matrix>
void qr_gs(Matrix & A, Matrix & Q, Matrix & R)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I m = A.shape[0];
    I n = A.shape[1];

    R[0][0] = A.col_norm_l2(0);

    for(I i = 0; i < m; i++)
        Q[i][0] = A[i][0]/R[0][0];

    for(I j = 1; j < n; j++)
    {
        for(I i = 0; i < m; i++)
            Q[i][j] = A[i][j];

        for(I k =0; k < j; k++)
        {
            for(I i = 0; i < m; i++)
            {
                R[k][j] += Q[i][k]*Q[i][j];
            }

            for(I i = 0; i < m; i++)
                Q[i][j] -= R[k][j]*Q[i][k]; 
        }

        R[j][j] = Q.col_norm_l2(j);

        for(I i = 0; i < m; i++)
            Q[i][j] /= R[j][j];
    }
}

template<typename Matrix>
void qr_gs(Matrix & A)
{
}

template<typename Matrix>
void qr_mgs(Matrix & A, Matrix & Q, Matrix & R)
{

}

template<typename Matrix>
void qr_mgs(Matrix & A)
{

}

} // end of AlgebraAlgrithom
} // end of WHYSC
#endif // end of linalg_h

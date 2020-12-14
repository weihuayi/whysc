#ifndef linalg_h
#define linalg_h

#include <iostream>
#include <initializer_list>

namespace WHYSC {
namespace AlgebraAlgrithom {

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

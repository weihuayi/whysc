#ifndef MatrixFactory_h
#define MatrixFactory_h

#include <iostream>
#include <cmath>
#include <initializer_list>
#include <vector>
#include <algorithm>

/*
 *
 * Notes:
 *  生成种测试用的矩阵和向量
 *
 */

namespace WHYSC {
namespace AlgebraObject {

template<typename Matrix>
void randn_matrix(Matrix & M, int low=0, int upper=10, int seed=0)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    I m = M.shape[0];
    I n = M.shape[1];

    if( seed == 0)
        srand((unsigned) time(NULL));
    else
        srand((unsigned) seed);

    for(I i = 0; i < m; i++)
    {
        for(I j = 0; j < n; j++)
        {
            int r = rand();
            M[i][j] = r%(upper - low) + low; 
        }
    }

    return;
}

template<typename I, typename SparseMatrix>
void laplace_1d(const I n, SparseMatrix & m)
{
    typedef typename SparseMatrix::Float  F;
    I nnz = 3*n-2;
    std::vector<F> data(nnz, 0.0);
    std::vector<I> row(nnz, 0);
    std::vector<I> col(nnz, 0);

    std::fill(data.begin(), data.begin() + n, 2.0);
    std::fill(data.begin() + n, data.end(), -1.0);


    // 对角线
    for(I i = 0; i < n; i++)
    {
        row[i] = i;
        col[i] = i;
    }

    // 副对角线
    for(I i = 0; i < n-1; i++)
    {
        row[n+i] = i;
        col[n+i] = i+1;

        row[2*n-1+i] = i+1;
        col[2*n-1+i] = i;
    }


    m.from_coo(n, n, nnz, row.data(), col.data(), data.data());
    return;
}

template<typename I, typename SparseMatrix>
void laplace_2d(const I N, SparseMatrix & m)
{
    typedef typename SparseMatrix::Float  F;
    I n = std::sqrt(N);
    I nnz = N + 4*(n-1)*n;
    std::vector<F> data(nnz, 0.0);
    std::vector<I> row(nnz, 0);
    std::vector<I> col(nnz, 0);

    std::fill(data.begin(), data.begin() + N, 4.0);
    std::fill(data.begin() + N, data.end(), -1.0);


    for(I i = 0; i < N; i++)
    {
        row[i] = i;
        col[i] = i;
    }

    I count = N;
    for(I i = 0; i < n; i++)
    {
        for(I j =0; j < n-1; j++)
        {
            row[count] = i*n + j;
            col[count] = row[count] + 1;
            
            col[count+1] = row[count];
            row[count+1] = col[count];
            count += 2;
        }
    }

    for(I j = 0; j < n; j++)
    {
        for(I i =0; i < n-1; i++)
        {
            row[count] = i*n + j;
            col[count] = row[count] + n;
            
            col[count+1] = row[count];
            row[count+1] = col[count];
            count += 2;
        }
    }

    m.from_coo(N, N, nnz, row.data(), col.data(), data.data());
    return;

}
/*
 *
 * Notes
 * -----
 *  M = n*n*n, 其中 n 为每个方向的节点个数
 *
 *  i: x 方向的指标
 *  j: y 方向的指标
 *  k: z 方向的指标
 *
 *  (i, j, k) 对应的节点的全局指标计算公式为
 *
 *      p = i*n*n + j*n + k
 */
template<typename I, typename SparseMatrix>
void laplace_3d(const I M, SparseMatrix & m)
{
    typedef typename SparseMatrix::Float  F;
    I n = std::pow(M, 1.0/3.0);
    I nnz = M + 6*n*n(n-1);
    std::vector<F> data(nnz, 0.0);
    std::vector<I> row(nnz, 0);
    std::vector<I> col(nnz, 0);

    std::fill(data.begin(), data.begin() + M, 6.0);
    std::fill(data.begin() + M, data.end(), -1.0);

    for(I i = 0; i < M; i++)
    {
        row[i] = i;
        col[i] = i;
    }

    I count = M;

    for(I i = 0; i < n; i++)
    {
        for( I j = 0; j < n; j++)
        {
            for(I k = 0; k < n-1; k++)
            {
                row[count] = i*n*n + j*n + k;
                col[count] = row[count] + 1;
                col[count+1] = row[count];
                row[count+1] = col[count];
                count += 2;
            }
        }
        for( I k = 0; k < n; k++)
        {
            for(I j = 0; j < n-1; j++)
            {
                row[count] = i*n*n + j*n + k;
                col[count] = row[count] + n;
                col[count+1] = row[count];
                row[count+1] = col[count];
                count += 2;
            }
        }
    }

    for(I j = 0; j < n; j++)
    {
        for( I k = 0; k < n; k++)
        {
            for(I i = 0; i < n-1; i++)
            {
                row[count] = i*n*n + j*n + k;
                col[count] = row[count] + n*n;
                col[count+1] = row[count];
                row[count+1] = col[count];
                count += 2;
            }
        }
    }
    m.from_coo(M, M, nnz, row.data(), col.data(), data.data());
    return;
}

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of MatrixFactory_h

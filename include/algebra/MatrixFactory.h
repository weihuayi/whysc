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

template<typename I, typename SparseMatrix>
void laplace_1d(const I n, SparseMatrix & m)
{
    typedef typename SparseMatrix::Float  F;
    std::vector<F> data(3*n-2, 0.0);
    std::vector<I> row(3*n-2, 0);
    std::vector<I> col(3*n-2, 0);

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


    m.from_coo(n, n, 3*n-2, row.data(), col.data(), data.data());
    return;
}

template<typename I, typename SparseMatrix>
void laplace_2d(const I N, SparseMatrix & m)
{
    typedef typename SparseMatrix::Float  F;
    I n = std::sqrt(N);
    I nnz = 4*(n-1)*n;
    std::vector<F> data(nnz, 0.0);
    std::vector<I> row(nnz, 0);
    std::vector<I> col(nnz, 0);

    std::fill(data.begin(), data.begin() + N, 4.0);
    std::fill(data.begin() + n, data.end(), -1.0);


    for(I i = 0; i < n; i++)
    {
        row[i] = i;
        col[i] = i;
    }

    for(I i = 0; i < n-1; i++)
    {
        row[n+i] = i;
        col[n+i] = i+1;

        row[2*n-1+i] = i+1;
        col[2*n-1+i] = i;
    }

    for(I i =0; i < 3*n-2; i++)
    {
        std::cout <<  row[i] << ", " << col[i] << ", " << data[i] << std::endl;
    }

    m.from_coo(n, n, 3*n-2, row.data(), col.data(), data.data());
    return;

}

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of MatrixFactory_h

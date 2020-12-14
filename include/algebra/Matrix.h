#ifndef Matrix_h
#define Matrix_h

#include <cmath>
#include <cassert>

#include <iostream>
#include <initializer_list>
#include <string>

namespace WHYSC {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct Matrix
{
    typedef F Float;
    typedef I Int;
    F ** data;
    I shape[2];

    static std::string format;

    /*
     * 默认构造函数
     */
    Matrix()
    {
        data = NULL;
        shape[0]=shape[1] = 0;
    }

    Matrix(I nr, I nc, F val=0.0)
    {
        shape[0] = nr;
        shape[1] = nc;
        data = new F*[shape[0]];
        for(I i=0; i < shape[0]; i++)
        {
            data[i] = new F[shape[1]];
            for(I j = 0; j < shape[1]; j++)
                data[i][j] = val;
        }
    }


    Matrix(const std::initializer_list<std::initializer_list<F>> &l)
    {
        shape[0] = l.size();
        shape[1] = l.begin()->size();

        data = new F*[shape[0]];
        I i = 0;
        for(auto & row: l)
        {
            data[i] = new F[shape[1]];
            I j=0;
            for(auto & val: row)
            {
                data[i][j] = val;
                j++;
            }
            i++;
        }
    }


    /*
     *
     * Notes
     * -----
     *  填充矩阵的对角线元素, 为一固定的值 val.
     *
     *  目前假设矩阵是方阵
     *
     */
    void fill_diag(const F val, const I diag=0)
    {
        I i=0;
        I j=0;

        if(diag >= 0)
        {
            j = diag;
            for(; j < shape[1]; i++, j++)
            {
                data[i][j] = val;
            }
        }
        else
        {
            i = std::abs(diag);
            for(; i < shape[0]; i++, j++)
            {
                data[i][j] = val;
            }
        }
    }

    void fill_diag(const F val[], const I diag=0)
    {
        I i=0;
        I j=0;

        if(diag >= 0)
        {
            j = diag;
            for(; j < shape[1]; i++, j++)
            {
                data[i][j] = val;
            }
        }
        else
        {
            i = std::abs(diag);
            for(; i < shape[0]; i++, j++)
            {
                data[i][j] = val;
            }
        }
    }

    inline F row_norm_l2(const I i)
    {
        F r = 0.0;
        for(I j = 0; j < shape[1]; j++)
            r += data[i][j]*data[i][j];
        r = std::sqrt(r);
        return r;
    }

    inline F col_norm_l2(const I j)
    {
        F c = 0.0;
        for(I i = 0; i < shape[0]; i++)
            c += data[i][j]*data[i][j];
        c = std::sqrt(c);
        return c;
    }

    ~Matrix()
    {
        if(data != NULL)
            delete [] data;
    }

    F & operator() (const I i, const I j)
    {
        return data[i][j];
    }

    const F & operator() (const I i, const I j) const
    {
        return data[i][j];
    }

    F * operator[](const I i) 
    {
        return data[i];
    }

    const F * operator[](const I i) const
    {
        return data[i];
    }

    F norm()
    {
        F sum = 0.0;
        for(auto i=0; i < shape[0]; i++)
            for(auto j=0; j < shape[1]; j++)
                sum += data[i][j]*data[i][j];
        return std::sqrt(sum);
    }

};

template<typename F, typename I>
std::string Matrix<F, I>::format = "full";

template<typename F, typename I>
inline Matrix<F, I> operator * (const Matrix<F, I> & m0, 
        const Matrix<F, I> & m1)
{
    auto nr = m0.shape[0];
    auto nc = m1.shape[1];
    auto n = m0.shape[1];
    Matrix<F, I> c(nr, nc);
    for(auto i=0; i < nr; i++)
    {
        for(auto j=0; j < nc; j++)
        {
            for(auto k=0; k < n; k++)
                c[i][j] += m0[i][k]*m1[k][j];
        }
    }

    return c;
}

template<typename F, typename I>
inline Matrix<F, I> operator - (const Matrix<F, I> & m0, 
        const Matrix<F, I> & m1)
{
    auto nr = m0.shape[0];
    auto nc = m0.shape[1];

    Matrix<F, I> c(nr, nc);
    for(auto i=0; i < nr; i++)
    {
        for(auto j=0; j < nc; j++)
        {
            c[i][j] = m0[i][j] - m1[i][j];
        }
    }

    return c;
}

template<typename F, typename I>
std::ostream& operator << (std::ostream & os, const Matrix<F, I> & m)
{
    std::cout << "Matrix("<< m.shape[0] << ","
        << m.shape[1] << "):" << std::endl;
    for(I i = 0; i < m.shape[0]; i ++)
    {
        for(I j = 0; j < m.shape[1]; j++)
        {
            os << m[i][j] << " ";
        }
        os << "\n";
    }
    os << "\n";
    return os;
}

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of Matrix_h

#ifndef MatrixFactory_h
#define MatrixFactory_h

#include <iostream>
#include <cmath>
#include <initializer_list>

/*
 *
 * Notes:
 *  生成种测试用的矩阵和向量
 *
 */

namespace WHYSC {
namespace AlgebraObject {

class MatrixFactory
{

template<class MatrixType>
static void laplace_1d(MatrixType & m)
{
    auto type = MatrixType::format;
    if(type == "full")
    {
    }
    else if(type == 'csr')
    {
    }
}

}; // end of class MatrixFactory

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of MatrixFactory_h

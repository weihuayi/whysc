#ifndef Vector_h
#define Vector_h

#include <algorithm>
#include <iostream>
#include <cmath>
#include <initializer_list>
#include <vector>

namespace WHYSC {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct Vector 
{
    typedef F Float;
    typedef I Int;

    F * data;
    I size;

    /*
     * 默认构造函数
     */
    Vector()
    {
        data = NULL;
        size = 0;
    }

    Vector(I s, F val=0.0)
    {
        size = s;
        init(val);
    }

    Vector(const std::initializer_list<F> &l)
    {
        size = l.size();
        data = new F[size];
        I i = 0;
        for(auto & val: l)
        {
            data[i] = val;
            i++;
        }
    }

    Vector<F, I>& operator = (const Vector<F, I> & rhs)
    {
        if(this != &rhs)
        {
            std::copy_n(rhs.data, size, data);
        }
        return *this;
    }

    void init(F val=0.0)
    {
        data = new F[size];
        for(I i=0; i < size; i++)
            data[i] = val;
    }

    ~Vector()
    {
        if(data != NULL)
            delete [] data;
    }

    F norm()
    {
        F sum = 0.0;
        for(I i=0; i < size; i++)
            sum += data[i]*data[i];
        return std::sqrt(sum);
    }

    F maxnorm()
    {
        F r = 0.0;
        for(I i=0; i < size; i++)
            if( r < std::abs(data[i]))
                r = std::abs(data[i]);
        return r;
    }

    F & operator[](const I i) 
    {
        return data[i];
    }

    const F & operator[](const I i) const
    {
        return data[i];
    }

};

template<typename F, typename I>
inline Vector<F, I> operator - (const Vector<F, I> & v0, 
        const Vector<F, I> & v1)
{
    Vector<F, I> r(v0.size);
    for(auto i=0; i < v0.size; i++)
    {
        r[i] = v0[i] - v1[i];
    }
    return r;
}

template<typename F, typename I>
std::ostream& operator << (std::ostream & os, const Vector<F, I> & v)
{
    std::cout << "Vector("<< v.size <<")" << std::endl;

    for(I i = 0; i < v.size; i++)
    {
        os << v[i] << " ";
    }
    os << std::endl;
    return os;
}

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of Vector_h

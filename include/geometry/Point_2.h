#ifndef Point_2_h
#define Point_2_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>
#include "Vector_2.h"

namespace WHYSC {

namespace GeometryObject {

template<typename F>
class Point_2 : public std::array<F, 2>
{
public:
    typedef typename std::array<F, 2> Base;
    typedef F Float;
    using Base::data;
public:

    Point_2()
    {
        std::fill_n(data(), 2, 0.0);
    }

    Point_2(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 2, data());
    }

    Point_2(F x, F y)
    {
        data()[0] = x;
        data()[1] = y;
    }

    template<typename P2> 
    Point_2(const P2 & p)
    {
        data()[0] = p[0];
        data()[1] = p[1];
    }

    Point_2(F * p)
    {
        data()[0] = p[0];
        data()[1] = p[1]; 
    }

    static int dimension() {return 2;}
    
    template<class V>
    Point_2<F> & operator += (const V & rhs)
    {
        for(int d = 0; d < 2; d++)
            data()[d] += rhs[d];
        return *this;
    }
    template<class V>
    Point_2<F> & operator -= (const V & rhs)
    {
        for(int d = 0; d < 2; d++)
            data()[d] -= rhs[d];
        return *this;
    }


};

template<typename F, typename V>
inline Point_2<F> operator + (const Point_2<F> & p, const V & v)
{
    Point_2<F> q;
    for(int d = 0; d < 2; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F>
inline Vector_2<F> operator - (const Point_2<F> & p, const Point_2<F> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Point_2<F> operator * (F w, const Point_2<F> & q)
{
    Point_2<F> p;
    for(int d = 0; d < 2; d++)
        p[d] = w*q[d];
    return p;
}

template<typename OS, typename F>
OS& operator << (OS & os, const Point_2<F> & p)
{
        return os << "Point_2(" << p[0] << ", " <<p[1] <<')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Point_2_h

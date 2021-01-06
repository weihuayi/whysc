#ifndef Point_3_h
#define Point_3_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>
#include "Vector_3.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Point_3 : public std::array<F, 3>
{
public:
    typedef typename std::array<F, 3> Base;
    typedef F Float;
    using Base::data;
public:

    Point_3()
    {
        data()[0] = 0.0;
        data()[1] = 0.0;
        data()[2] = 0.0;
    }

    Point_3(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 3, data());
    }

    Point_3(F x, F y, F z)
    {
        data()[0] = x;
        data()[1] = y;
        data()[2] = z;
    }

    template<typename P3>
    Point_3(const P3 & p)
    {
        data()[0] = p[0];
        data()[1] = p[1];
        data()[2] = p[2];
    }

    Point_3(F * p)
    {
        data()[0] = p[0];
        data()[1] = p[1]; 
        data()[2] = p[2]; 
    }

    static int dimension() {return 3;}

    template<class V>
    Point_3<F> & operator -= (const V & rhs)
    {
        for(int d = 0; d < 3; d++)
            data()[d] -= rhs[d];
        return *this;
    }

    template<class V>
    Point_3<F> & operator += (const V & rhs)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] += rhs[d];
        return *this;
    }
};

template<typename F, typename V>
inline Point_3<F> operator + (const Point_3<F> & p, const V & v)
{
    Point_3<F> q;
    for(int d = 0; d < 3; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F>
inline Vector_3<F> operator - (const Point_3<F> & p, const Point_3<F> & q)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Point_3<F> operator * (F w, const Point_3<F> & q)
{
    Point_3<F> p;
    for(int d = 0; d < 3; d++)
        p[d] = w*q[d];
    return p;
}

template<typename OS, typename F>
OS& operator << (OS & os, const Point_3<F> & p)
{
    return os << "Point_3(" << p[0] << ", " << p[1] << ", " << p[2] << ')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Point_3_h

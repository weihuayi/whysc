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
class Point_2 
{
public:
    typedef F Float;
public:

    Point_2()
    {
       m_data[0] = 0.0;
       m_data[1] = 0.0;
    }

    Point_2(const std::initializer_list<F> &l)
    { 
        m_data[0] = l.begin()[0];
        m_data[1] = l.begin()[1];
    }

    Point_2(F x, F y)
    {
        m_data[0] = x;
        m_data[1] = y;
    }

    template<typename P2> 
    Point_2(const P2 & p)
    {
        m_data[0] = p[0];
        m_data[1] = p[1];
    }

    Point_2(F * p)
    {
        m_data[0] = p[0];
        m_data[1] = p[1]; 
    }

    Point_2& operator=(const Point_2& p)
    {
        m_data[0] = p[0];
        m_data[1] = p[1];
        return *this;
    }

    static int size() {return 2;}

    static int dimension() {return 2;}

    F & x()
    {
        return m_data[0];
    }

    F & y()
    {
        return m_data[1];
    }
    
    F & operator [] (const int i)
    {
        return m_data[i];
    }

    const F & operator [] (const int i) const
    {
        return m_data[i];
    }

    template<class V>
    Point_2<F> & operator += (const V & rhs)
    {
        m_data[0] += rhs[0];
        m_data[1] += rhs[1]; 
        return *this;
    }
    template<class V>
    Point_2<F> & operator -= (const V & rhs)
    {
        for(int d = 0; d < 2; d++)
            m_data[d] -= rhs[d];
        return *this;
    }

private:
    F m_data[2];
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

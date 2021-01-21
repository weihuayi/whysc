#ifndef Point_3_h
#define Point_3_h

#include <initializer_list>
#include <assert.h>
#include "Vector_3.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Point_3 
{
public:
    typedef F Float;
public:

    Point_3()
    {
        m_data[0] = 0.0;
        m_data[1] = 0.0;
        m_data[2] = 0.0;
    }

    Point_3(const std::initializer_list<F> &l)
    { 
        m_data[0] = l.begin()[0];
        m_data[1] = l.begin()[1];
        m_data[2] = l.begin()[2];
    }

    Point_3(F x, F y, F z)
    {
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
    }

    template<typename P3>
    Point_3(const P3 & p)
    {
        m_data[0] = p[0];
        m_data[1] = p[1];
        m_data[2] = p[2];
    }

    Point_3(F * p)
    {
        m_data[0] = p[0];
        m_data[1] = p[1]; 
        m_data[2] = p[2]; 
    }

    static int size() {return 3;}
    static int dimension() {return 3;}

    F & x()
    {
        return m_data[0];
    }

    F & y()
    {
        return m_data[1];
    }

    F & z()
    {
        return m_data[2];
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
    Point_3<F> & operator -= (const V & rhs)
    {
        m_data[0] -= rhs[0];
        m_data[1] -= rhs[1];
        m_data[2] -= rhs[2];
        return *this;
    }

    template<class V>
    Point_3<F> & operator += (const V & rhs)
    {
        m_data[0] += rhs[0];
        m_data[1] += rhs[1];
        m_data[2] += rhs[2];
        return *this;
    }


private:
    F m_data[3];
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

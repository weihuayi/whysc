#ifndef Vector_3_h
#define Vector_3_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Vector_3
{
public:
    typedef F Float;
public:

    Vector_3()
    {
       m_data[0] = 0.0;
       m_data[1] = 0.0;
       m_data[2] = 0.0;
    }

    Vector_3(const std::initializer_list<F> &l)
    { 
        for(int d = 0; d < 3; d++)
           m_data[d] = l.begin()[d];
    }

    Vector_3(const Vector_3 &v)
    { 
        for(int d = 0; d < 3; d++)
           m_data[d] = v[d];
    }

    Vector_3(F vx, F vy, F vz)
    {
       m_data[0] = vx;
       m_data[1] = vy;
       m_data[2] = vz;
    }

    static int size() {return 3;}
    static int dimension() {return 3;}

    F & operator [] (const int i)
    {
        return m_data[i];
    }

    const F & operator [] (const int i) const
    {
        return m_data[i];
    }

    F squared_length()
    {

        return m_data[0]*m_data[0] + m_data[1]*m_data[1] +  m_data[2]*m_data[2];
    }

    template<class RVector_3>
    F operator * (const RVector_3 & w)
    {
        return m_data[0]*w[0] + m_data[1]*w[1] +  m_data[2]*w[2];
    }

    template<class RVector_3>
    F dot (const RVector_3 & w)
    {
        return m_data[0]*w[0] + m_data[1]*w[1] +  m_data[2]*w[2];
    }

    Vector_3<F> & operator *= (const F & s)
    {
        for(int d = 0; d < 3; d++)
            m_data[d] *= s;
        return *this;
    }


    Vector_3<F> & operator /= (const F & s)
    {
        for(int d = 0; d < 3; d++)
            m_data[d] /= s;
        return *this;
    }

    template<class RVector_3>
    Vector_3<F> & operator += (const RVector_3 & w)
    {
        for(int d = 0; d < 3; d++)
            m_data[d] += w[d];
        return *this;
    }

    template<class RVector_3>
    Vector_3<F> & operator -= (const RVector_3 & w)
    {
        for(int d = 0; d < 3; d++)
            m_data[d] -= w[d];
        return * this;
    }
private:
    F m_data[3];
};

template<typename F>
inline Vector_3<F> operator + (const Vector_3<F> & v0, const Vector_3<F> & v1)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = v0[d] + v1[d];
    return v;
}

template<typename F>
inline Vector_3<F>  cross(const Vector_3<F> & v, const Vector_3<F> & w)
{
    Vector_3<F> r;
    r[0] = v[1]*w[2] - v[2]*w[1];
    r[1] = v[2]*w[0] - v[0]*w[2];
    r[2] = v[0]*w[1] - v[1]*w[0];
    return r;
}

template<typename F>
inline F  dot(const Vector_3<F> & v, const Vector_3<F> & w)
{
    return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

template<typename F>
inline Vector_3<F> operator - (const Vector_3<F> & v0)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = -v0[d];
    return v;
}


template<typename F>
inline Vector_3<F> operator - (const Vector_3<F> & v0, const Vector_3<F> & v1)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = v0[d] - v1[d];
    return v;
}

template<typename F>
inline Vector_3<F> operator * (const F w, const Vector_3<F> & v0)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = w*v0[d];
    return v;
}

template<typename F>
inline Vector_3<F> operator * (const Vector_3<F> & v0, const F w)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = w*v0[d];
    return v;
}

template<typename F, typename I>
inline Vector_3<F> operator / (const Vector_3<F> & v0, const I w)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = v0[d]/w;
    return v;
}

template<typename OS, typename F>
OS& operator << (OS & os, const Vector_3<F> & v)
{
        return os << "Vector_3(" << v[0] << ", " << v[1] << ", " << v[2] << ')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Vector_3_h

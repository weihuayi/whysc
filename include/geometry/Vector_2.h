#ifndef Vector_2_h
#define Vector_2_h

#include <initializer_list>
#include <assert.h>

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Vector_2
{
public:
    typedef F Float;
public:

    Vector_2()
    {
       m_data[0] = 0.0;
       m_data[1] = 0.0;
    }

    Vector_2(const std::initializer_list<F> &l)
    { 
        m_data[0] = l.begin()[0];
        m_data[1] = l.begin()[1];
    }

    Vector_2(const Vector_2 &v)
    { 
        m_data[0] = v[0];
        m_data[1] = v[1];
    }

    Vector_2(F vx, F vy)
    {
        m_data[0] = vx;
        m_data[1] = vy;
    }

    Vector_2& operator=(const Vector_2& v)
    {
        m_data[0] = v[0];
        m_data[1] = v[1];
        return *this;
    }

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

    static int size() {return 2;}
    static int dimension() {return 2;}

    F squared_length()
    {
        return m_data[0]*m_data[0] + m_data[1]*m_data[1];
    }

    template<class RVector_2>
    F dot(const RVector_2 & w)
    {
        return m_data[0]*w[0] + m_data[1]*w[1];
    }

    template<class RVector_2>
    F cross(const RVector_2 & w)
    {
        return m_data[0]*w[1] - w[0]*m_data[1];
    }

    template<class RVector_2>
    F operator * (const RVector_2 & w)
    {
        return m_data[0]*w[0] + m_data[1]*w[1];
    }

    Vector_2<F> & operator *= (const F & s)
    {
        m_data[0] *= s;
        m_data[1] *= s;
        return *this;
    }


    Vector_2<F> & operator /= (const F & s)
    {
        m_data[0] /= s;
        m_data[1] /= s;
        return *this;
    }

    template<class RVector_2>
    Vector_2<F> & operator += (const RVector_2 & w)
    {
        m_data[0] += w[0];
        m_data[1] += w[1];
        return * this;
    }


    template<class RVector_2>
    Vector_2<F> & operator -= (const RVector_2 & w)
    {
        m_data[0] -= w[0];
        m_data[1] -= w[1];
        return * this;
    }

private:
    F m_data[2];
};

template<typename F>
inline F  dot(const Vector_2<F> & v, const Vector_2<F> & w)
{
    return v[0]*w[0] + v[1]*w[1];
}

template<typename F>
inline F cross(const Vector_2<F> & v0, const Vector_2<F> & v1)
{
    return v0[0]*v1[1] - v1[0]*v0[1];
}

template<typename F>
inline Vector_2<F> operator + (const Vector_2<F> & v0, const Vector_2<F> & v1)
{
    Vector_2<F> v;
    for(auto d = 0; d < 2; d++)
        v[d] = v0[d] + v1[d];
    return v;
}

template<typename F>
inline Vector_2<F> operator - (const Vector_2<F> & v0, const Vector_2<F> & v1)
{
    Vector_2<F> v;
    for(auto d = 0; d < 2; d++)
        v[d] = v0[d] - v1[d];
    return v;
}

template<typename F>
inline Vector_2<F> operator * (F w, const Vector_2<F> & v0)
{
    Vector_2<F> v;
    for(auto d = 0; d < 2; d++)
        v[d] = w*v0[d];
    return v;
}

template<typename F>
inline Vector_2<F> operator / (const Vector_2<F> & v0, const F w)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = v0[d]/w;
    return v;
}

template<typename F>
inline Vector_2<F> operator - (const Vector_2<F> & v0)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = -v0[d];
    return v;
}

template<typename OS, typename F>
OS& operator << (OS & os, const Vector_2<F> & v)
{
        return os << "Vector_2(" << v[0] << ", " << v[1] <<')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Vector_2_h

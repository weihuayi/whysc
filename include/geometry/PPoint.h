#ifndef PPoint_h
#define PPoint_h

#include "Vector_2.h"
#include "Vector_3.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F, int DIM>
class PPoint
{
public:
    typedef F Float;
private:
    F * m_data; // 注意这里 PPoint 类不负责数据的管理
public:

    PPoint()
    {
        m_data = NULL;
    }

    PPoint(F * p)
    {
        m_data = p;
    }

    PPoint & operator=(F * p)
    {
        this->m_data = p;
        return *this;
    }

    PPoint(const PPoint & p)
    {
        m_data = p.m_data;
    }
    PPoint& operator=(const PPoint& p)
    {
        m_data = p.mdata;
    }

    static int dimension() {return DIM;}

    F & operator [] (const int i)
    {
        return m_data[i];
    }

    const F & operator [] (const int i) const
    {
        return m_data[i];
    }

    template<class V>
    PPoint<F, DIM> & operator += (const V & rhs)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] += rhs[d];
        return *this;
    }

    template<class V>
    PPoint<F, DIM> & operator -= (const V & rhs)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] -= rhs[d];
        return *this;
    }

    PPoint<F, DIM> & operator *= (const F w)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] *= w;
        return *this;
    }

    PPoint<F, DIM> & operator /= (const F w)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] /= w;
        return *this;
    }

};

template<typename F>
inline Vector_2<F> operator - (const PPoint<F, 2> & p, const PPoint<F, 2> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Vector_3<F> operator - (const PPoint<F, 3> & p, const PPoint<F, 3> & q)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename OS, typename F, int DIM>
OS& operator << (OS & os, const PPoint<F, DIM> & p)
{
    if(DIM == 2)
        return os << "PPoint_2(" << p[0] << ", " <<p[1] <<')';
    else if(DIM == 3)
        return os << "PPoint_3(" << p[0] << ", " <<p[1] <<", " << p[2]<< ')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of PPoint_h

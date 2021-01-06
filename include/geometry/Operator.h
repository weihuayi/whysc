#ifndef Operator_h
#define Operator_h

#include "Vector_2.h"
#include "Vector_3.h"
#include "Point_2.h"
#include "Point_3.h"
#include "PPoint.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F>
inline Point_2<F> operator * (F w, const PPoint<F, 2> & q)
{
    Point_2<F> p;
    for(int d = 0; d < 2; d++)
        p[d] = w*q[d];
    return p;
}

template<typename F>
inline Point_3<F> operator * (F w, const PPoint<F, 3> & q)
{
    Point_3<F> p;
    for(int d = 0; d < 3; d++)
        p[d] = w*q[d];
    return p;
}

template<typename F>
inline Point_2<F> operator + (const PPoint<F, 2> & p, const Vector_2<F> & v)
{
    Point_2<F> q;
    for(int d = 0; d < 2; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F>
inline Vector_2<F> operator - (const PPoint<F, 2> & p, const Point_2<F> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Vector_2<F> operator - (const Point_2<F> & p, const PPoint<F, 2> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Vector_3<F> operator - (const PPoint<F, 3> & p, const Point_3<F> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Vector_3<F> operator - (const Point_3<F> & p, const PPoint<F, 3> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = p[d] - q[d];
    return v;
}

} // end of namespace GeometryObject
} // end of namespace WHYSC

#endif // end of Operator_h

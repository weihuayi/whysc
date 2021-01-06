#ifndef Bisection_alg_h
#define Bisection_alg_h

#include <cmath>
namespace WHYSC {

namespace GeoAlg {

template<class GK> 
class Bisection_alg
{
public:
    typedef typename GK::Float Float;
    typedef typename GK::Int Int;
public:
    Bisection_alg(){}

    template<class LevelSetFunction, class Point >
    Point operator () (LevelSetFunction & fun, const Point & p1, const Point & p2) 
    {
        Point a = p1;
        Point b = p2;
        Int a_sign = GK::sign(fun(a));
        Int b_sign = GK::sign(fun(b));

        Float h = std::sqrt((b - a).squared_length());
        Point m = GK::midpoint(a,b); 
        Int m_sign = GK::sign(fun(m));
        while(h > GK::eps())
        {
            if(m_sign == 0)
            {
                return m; 
            }
            else if( a_sign*m_sign < 0)
            {
                b = m;
                b_sign = m_sign;
            }
            else if( b_sign*m_sign < 0)
            {
                a = m;
                a_sign = m_sign;
            }

            h = std::sqrt((b - a).squared_length()); 

            m = GK::midpoint(a,b); 
            m_sign = GK::sign(fun(m));
        }

        return m;
    }
};

} // end of namespace GeoAlg

} // end of namespace WHYSC
#endif // end of Bisection_alg_h

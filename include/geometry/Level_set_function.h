#ifndef Level_set_function_h
#define Level_set_function_h

#include <vector>
#include <algorithm>
#include <cmath>

namespace WHYSC {

namespace LevelSetFunction {


template<class GK>
class Circle_2
{
public:
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Float Float;
    typedef typename GK::Int Int;
public:
    /*
     * Constructor
     *
     */
    Circle_2(): _center(0, 0), _r(1)
    {
    }

    Circle_2(const Float x, const Float y, const Float r): _center(x, y), _r(r)
    {
    }

    Circle_2(const Circle_2 & c)
    {
        _center = c._center;
        _r = c._r;
    }

    Float operator () (const Point_2 & p)
    {

        return  (p[0]-_center[0])*(p[0]-_center[0]) + (p[1]-_center[1])*(p[1]-_center[1]) - _r*_r;
    }


    Int sign(const Point_2 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Point_2 center() { return _center;}
    Float radius() { return _r;}

private:
    Point_2 _center;
    Float _r;
    
};

template<class GK>
class Signed_distance_circle_2
{
public:
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Float Float;
    typedef typename GK::Int Int;
public:
    /*
     * Constructor
     *
     */
    Signed_distance_circle_2(): _center(Float(0), Float(0)), _r(Float(1))
    {
    }

    Signed_distance_circle_2(const Float x, const Float y, const Float r): _center(x, y), _r(r)
    {
    }

    Signed_distance_circle_2(const Signed_distance_circle_2 & c)
    {
        _center = c._center;
        _r = c._r;
    }

    Float operator () (const Point_2 & p)
    {

        return  std::sqrt((p-_center).squared_length()) - _r;
    }

    Int sign(const Point_2 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Point_2 point(const double * p)
    {
        return GK::point_2(p);
    }

private:
    Point_2 _center;
    Float _r;
    
};


template<class GK>
class Sphere_3
{
public:
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Float Float;
    typedef typename GK::Int   Int;
public:
    Sphere_3():_center(Float(0.0), Float(0.0), Float(0.0)), _r(Float(1)) {}
    Sphere_3(const Float x, const Float y, const Float z, const Float r):_center(x, y, z), _r(r){}
    Sphere_3(const Sphere_3 & s)
    {
        _center = s._center;
        _r = s._r;
    }

    Float operator () (const Point_3 & p)
    {

        return  (p-_center).squared_length() - _r*_r;
    }

    Int sign(const Point_3 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Point_3 point(const Float *p)
    {
        return GK::point_3(p);
    }
    Point_3 center() { return _center;}
    Float radius() { return _r;}

private:
    Point_3 _center;
    Float _r;
};

template<class GK>
class Double_torus_3 
{
public:
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Float Float;
    typedef typename GK::Int   Int;
public:

    Float operator () (const Point_3 & p)
    {

        auto x2 = p[0]*p[0];
        auto y2 = p[1]*p[1];
        auto y4 = y2*y2;
        auto z2 = p[2]*p[2];
        
        return  x2*(x2 - 1)*(x2*(x2-1) + 2*y2) + y4 + z2 - 0.04;
    }

    Int sign(const Point_3 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Point_3 point(const Float *p)
    {
        return GK::point_3(p);
    }
};

template<class GK>
class Orthocircle_3 
{
public:
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Float Float;
    typedef typename GK::Int   Int;
public:

    Float operator () (const Point_3 & p)
    {

        auto x2 = p[0]*p[0];
        auto y2 = p[1]*p[1];
        auto y4 = y2*y2;
        auto z2 = p[2]*p[2];

        auto r0 = x2 + y2 + z2;
        auto r1 = x2 + y2 - 1;
        auto r2 = y2 + z2 - 1;
        auto r3 = z2 + x2 - 1;
        
        return  (r1*r1 + z2)*(r2*r2 + x2)*(r3*r3 + y2) - 0.075*0.075*(1 + 3*r0);
    }

    Int sign(const Point_3 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Point_3 point(const Float *p)
    {
        return GK::point_3(p);
    }
};

template<class GK>
class Signed_distance_sphere_3
{
public:
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Vector_3 Vector_3;
    typedef typename GK::Float Float;
    typedef typename GK::Int   Int;
public:
    Signed_distance_sphere_3():_center(Float(0.0), Float(0.0), Float(0.0)), _r(Float(1.0)) {}
    Signed_distance_sphere_3(const Float x, const Float y, const Float z, const Float r):_center(x, y, z), _r(r){}
    Signed_distance_sphere_3(const Point_3 & c, const Float r):_center(c), _r(r){}
    Signed_distance_sphere_3(const Signed_distance_sphere_3 & s)
    {
        _center = s._center;
        _r = s._r;
    }

    Float operator () (const Point_3 & p)
    {
        return  std::sqrt((p-_center).squared_length()) - _r;
    }

    Int sign(const Point_3 & p)
    {
        Float val = this->operator () (p);
        return GK::sign(val);
    }

    Vector_3 gradient(const Point_3 & p)
    {
        Float d = std::sqrt((p - _center).squared_length());
        Vector_3 v = (p - _center)/d;
        return v;
    }

    Point_3 project(const Point_3 & p)
    {
        Float d = this->operator ()(p);
        Vector_3 v = gradient(p);
        Point_3 q = p - d*v;
        return q;
    }

    Point_3 center() { return _center;}
    Float radius() { return _r;}
private:
    Point_3 _center;
    Float _r;
};


//template<class GK>
//class Union
//{
//public:
//    Union(Base & lsf_1, Base & lsf_2)
//    {
//        _lsfs.push_back(lsf_1);
//        _lsfs.push_back(lsf_2);
//    }
//
//    Union(const std::vector<Base> & lsfs): _lsfs(lsfs)
//    {
//
//    }
//
//    double operator () (const Point & p)
//    {
//        int n = _lsfs.size();
//        std::vector<double> vals(n);
//        for(int i = 0; i < n; i++)
//            vals[i] = _lsfs[i](p);
//
//        double max_val = std::min_element(vals.begin(), vals.end());
//        return max_val;
//    }
//
//    int sign(const Point & p)
//    {
//        double val = this->operator () (p);
//        return Base::sign(val);
//    }
//private:
//    std::vector<Base> _lsfs;
//};

} // end of namespace LevelSetFunction

} // end of namespace WHYSC
#endif // end of Level_set_function_h

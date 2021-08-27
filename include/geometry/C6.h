#ifndef C6_h
#define C6_h

#include <map>
#include <vector>
#include <math.h>
#include <iostream>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class C6
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;

public:
  C6(): m_r(1.0)
  {
    double a = std::sqrt(3);
    double b = std::sqrt(3)/2.0;
    m_center.resize(6);
    m_center[0] = Point(-a, 0.0, 0.0);
    m_center[1] = Point(-b, -1.5, 0.0);
    m_center[2] = Point(b, -1.5, 0.0);
    m_center[3] = Point(a, 0.0, 0.0);
    m_center[4] = Point(b, 1.5, 0.0);
    m_center[5] = Point(-b, 1.5, 0.0);

    m_norm.resize(6);
    m_norm[0] = Vector(-0.5, -b, 0.0);
    m_norm[1] = Vector(0.5, -b, 0.0);
    m_norm[2] = Vector(1.0, 0.0, 0.0);
    m_norm[3] = Vector(0.5, b, 0.0);
    m_norm[4] = Vector(-0.5, b, 0.0);
    m_norm[5] = Vector(-1.0, 0.0, 0.0);
  }
                                           
  void project_to_face(const int fid, Point & p)
  {
    if(fid<7)
    {
      auto center = m_center[fid-1];
      auto v = p - center;
      v = (m_r/std::sqrt(v.squared_length()))*v;
      p = center+v;
    }
    else
    {
      Point O(0, 0, 0);
      auto n = m_norm[fid-7];
      auto v = p-O;
      v = v - dot(v, n)*n;
      p = O+v;
    }
  }

  void project_to_edge(const int eid, Point & p) {}

  void project_vector_to_face(const int fid, const Point p, Vector & v)
  {
    if(fid<7)
    {
      auto center = m_center[fid-1];
      auto n = p - center;
      v = v - dot(v, n)*n/n.squared_length();
    }
    else
    {
      Vector O(0, 0, 0);
      auto n = m_norm[fid-7];
      v = v - dot(v, n)*n;
    }
  }

  void project_vector_to_edge(const int eid, const Point p, Vector & v) {}

private:
  double m_r;
  std::vector<Point> m_center;
  std::vector<Vector> m_norm;
};

}
}


#endif // end of C6_h

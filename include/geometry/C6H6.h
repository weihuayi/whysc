#ifndef C6H6_h
#define C6H6_h

#include <map>
#include <vector>
#include <math.h>
#include <iostream>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class C6H6
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;

public:
  C6H6(): m_r(1.1*std::sqrt(2))
  {
    double g3 = std::sqrt(3)/2;
    m_center.resize(6);
    m_center[0] = 4.0*Point(1.0, 0.0, 0.0);
    m_center[1] = 4.0*Point(0.5, g3, 0.0);
    m_center[2] = 4.0*Point(-0.5, g3, 0.0);
    m_center[3] = 4.0*Point(-1.0, 0.0, 0.0);
    m_center[4] = 4.0*Point(-0.5, -g3, 0.0);
    m_center[5] = 4.0*Point(0.5, -g3, 0.0);
  }
                                           
  void project_to_face(const int fid, Point & p)
  {
    auto center = m_center[fid];
    auto v = p - center;
    v = (m_r/std::sqrt(v.squared_length()))*v;
    p = center+v;
  }

  void project_to_edge(const int eid, Point & p) {}

  void project_vector_to_face(const int fid, const Point p, Vector & v)
  {
    auto center = m_center[fid];
    auto n = p - center;
    v = v - dot(v, n)*n/n.squared_length();
  }

  void project_vector_to_edge(const int eid, const Point p, Vector & v) {}

private:
  double m_r;
  std::vector<Point> m_center;
};

}
}


#endif // end of C6H6_h

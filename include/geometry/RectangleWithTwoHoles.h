#ifndef RectangleWithTwoHoles_h
#define RectangleWithTwoHoles_h

#include <math.h>
#include <vector>
#include <map>

namespace WHYSC {
namespace GeometryModel{

template<typename GK>
class RectangleWithTwoHoles
{
public:
  typedef typename GK::Point_2 Point;
  typedef typename GK::Vector_2 Vector;
  typedef typename std::array<int, 2> Line;

public:
  RectangleWithTwoHoles(): m_r(2),  m_center(2), m_point(4), m_line(4)
  {
    m_r[0] = 0.5;
    m_r[1] = 1.0;

    double g2 = std::sqrt(2)/2.0;
    m_center[0] = Point(0.8, 0.8);
    m_center[1] = Point(1.6+g2, 1.6+g2);

    m_point[0] = Point(0.0, 0.0);
    m_point[1] = Point(4.0, 0.0);
    m_point[2] = Point(4.0, 4.0);
    m_point[3] = Point(0.0, 4.0);

    m_line[0] = Line({0, 1});
    m_line[1] = Line({2, 3});
    m_line[2] = Line({1, 2});
    m_line[3] = Line({3, 0});
  }

  template<typename T>
  void project_to_edge(const int & eid, T & p)
  {
    if(eid<3)
    {
      const auto & center = m_center[eid-1];
      const auto & r = m_r[eid-1];
      auto v = p - center;
      v = r*(v/std::sqrt(v.squared_length()));
      p = center+v;
    }
    else
    {
      const auto & p0 = m_point[m_line[eid-3][0]];
      const auto & p1 = m_point[m_line[eid-3][1]];
      auto v = p1 - p0;
      auto v0 = p - p0;

      double k = dot(v0, v)/dot(v, v);
      p[0] = p0[0] + k*v[0];
      p[1] = p0[1] + k*v[1];
    }
  }

  void project_to_face(const int eid, Point & p) {}

  void project_vector_to_face(const int fid, const Point p, Vector & v) {}

  void project_vector_to_edge(const int eid, const Point p, Vector & v) 
  {
    if(eid<3)
    {
      const auto & center = m_center[eid-1];
      const auto & r = m_r[eid-1];
      auto d = p - center;
      v = v - dot(d, v)*d/dot(d, d);
    }
    else
    {
      const auto & p0 = m_point[m_line[eid-3][0]];
      const auto & p1 = m_point[m_line[eid-3][1]];
      auto d = p1 - p0;
      v = dot(d, v)*d/dot(d, d);
    }
  }

private:
  std::vector<double> m_r;
  std::vector<Point> m_center;
  std::vector<Line> m_line;
  std::vector<Point> m_point;
};

} // end of GeometryModel

} // end of WHYSC

#endif // end of RectangleWithTwoHoles_h

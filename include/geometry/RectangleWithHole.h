#ifndef RectangleWithHole_h
#define RectangleWithHole_h

#include <math.h>
#include <vector>
#include <map>

namespace WHYSC {
namespace GeometryModel{

template<typename GK>
class RectangleWithHole
{
public:
  typedef typename GK::Point_3 Point;

  typedef typename std::vector<int> Line;
  typedef typename std::vector<int> Face;
  typedef typename std::vector<int> Volume;
  class Circle
  {
  public:
    Circle(){}
    double radius;
    Point center;
  };
public:
  RectangleWithHole()
  {
    add_point({0, 0, 0}, 1);
    add_point({1, 0, 0}, 2);
    add_point({1, 1, 0}, 3);
    add_point({0, 1, 0}, 4);

    add_line({1, 2}, 1);
    add_line({2, 3}, 2);
    add_line({3, 4}, 3);
    add_line({4, 1}, 4);

    add_circle({0.5, 0.5, 0}, 0.3, 1);
    add_face({1, 2, 3, 4, -5, -6, -7, -8}, 1);
  }
  void add_point(std::initializer_list<double> point, int tag)
  {
    m_points[tag] = Point(point); 
  }

  void add_line(std::initializer_list<int> line, int tag)
  {
    m_lines[tag] = Line(line);
  }

  void add_face(std::initializer_list<int> face, int tag)
  {
    m_faces[tag] = Face(face);
  }

  void add_circle(std::initializer_list<double> point, double r, int tag)
  {
    Circle c;
    c.center = Point(point);
    c.radius = r;
    m_circles[tag] = c; 
  }

  void project_to_edge(const int eid, Point & p)
  {
    if(eid<5)
    {
      auto & p0 = m_points[m_lines[eid][0]];
      auto & p1 = m_points[m_lines[eid][1]];
      auto v = p1 - p0;

      double k = dot(p-p0, v)/dot(v, v);
      p = p0 + k*v;
    }
    else
    {
      int cirid = (eid-5)/4;
      auto center = m_circles[cirid].center;
      auto r = m_circles[cirid].radius;
      auto v = p - center;
      v = r*v/std::sqrt(v.squared_length());
      p = v+center;
    }
  }

  std::map<int, Point> & get_points()
  {
    return m_points;
  }

  std::map<int, Line> & get_lines()
  {
    return m_lines;
  }

  std::map<int, Face> & get_faces()
  {
    return m_faces;
  }
 
  std::map<int, Circle> & get_circles()
  {
    return m_circles;
  }

private:
  std::map<int, Point> m_points;
  std::map<int, Line> m_lines;
  std::map<int, Face> m_faces;
  std::map<int, Volume> m_volumes;
  std::map<int, Circle> m_circles;
};

} // end of GeometryModel

} // end of WHYSC

#endif // end of RectangleWithHole_h

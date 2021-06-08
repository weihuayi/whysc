#ifndef CubeModel_h
#define CubeModel_h

#include <map>
#include <vector>
#include <math.h>
#include <iostream>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class CubeModel
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;

  typedef typename std::vector<int> Line;
  typedef typename std::vector<int> Face;
  typedef typename std::vector<int> Cube;
  typedef typename GK::Sphere_3 Sphere;

public:
  CubeModel()
  {
    add_point({0.0, 0.0, 0.0}, 1);
    add_point({1.0, 0.0, 0.0}, 2);
    add_point({1.0, 1.0, 0.0}, 3);
    add_point({0.0, 1.0, 0.0}, 4);
    add_point({0.0, 0.0, 1.0}, 5);
    add_point({1.0, 0.0, 1.0}, 6);
    add_point({1.0, 1.0, 1.0}, 7);
    add_point({0.0, 1.0, 1.0}, 8);
                                           
    add_line({1, 2}, 1);
    add_line({2, 3}, 2);
    add_line({3, 4}, 3);
    add_line({4, 1}, 4);
    add_line({5, 6}, 5);
    add_line({6, 7}, 6);
    add_line({7, 8}, 7);
    add_line({8, 5}, 8);
    add_line({1, 5}, 9);
    add_line({2, 6}, 10);
    add_line({3, 7}, 11);
    add_line({4, 8}, 12);
                                           
    add_face({-1, -4, -3, -2}, 1);
    add_face({5, 6, 7, 8}, 2);
    add_face({1, 10, -5, -9}, 3);
    add_face({3, 12, -7, -11}, 4);
    add_face({2, 11, -6, -10}, 5);
    add_face({4, 9, -8, -12}, 6);
                                           
    add_volume({1, 2, 3, 4, 5, 6}, 1);
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

  void add_volume(std::initializer_list<int> vol, int tag)
  {
    m_volumes[tag] = Cube(vol);
  }

  bool sphere_flag()
  {
    return false;
  }

  std::map<int, Sphere> & get_spheres()
  {
    return m_spheres;
  }

  void project_to_face(const int fid, Point & p)
  {
    Point p0, p1, p2;
    if(m_faces[fid][0] > 0)
      p0 = m_points[m_lines[m_faces[fid][0]][0]];
    else
      p0 = m_points[m_lines[-m_faces[fid][0]][1]];

    if(m_faces[fid][1] > 0)
      p1 = m_points[m_lines[m_faces[fid][1]][0]];
    else
      p1 = m_points[m_lines[-m_faces[fid][1]][1]];

    if(m_faces[fid][2] > 0)
      p2 = m_points[m_lines[m_faces[fid][2]][0]];
    else
      p2 = m_points[m_lines[-m_faces[fid][2]][1]];

    auto v0 = p1 - p0;
    auto v1 = p2 - p1;
    auto v2 = p - p0;
    double a = dot(v0, v0);
    double b = dot(v0, v1);
    double c = dot(v1, v1);
    double f = dot(v2, v0);
    double g = dot(v2, v1);

    double k = (f*c-b*g)/(a*c-b*b);
    double m = (a*g-b*f)/(a*c-b*b);
    p = p0 + k*v0 + m*v1;
  }

  void project_to_edge(const int eid, Point & p)
  {
    auto & p0 = m_points[m_lines[eid][0]];
    auto & p1 = m_points[m_lines[eid][1]];
    auto v = p1 - p0;

    double k = dot(p-p0, v)/dot(v, v);
    p = p0 + k*v;
  }

  void get_point_normal(const int fid, const Point p, Vector & n)
  {
    Point p0, p1, p2;
    if(m_faces[fid][0] > 0)
      p0 = m_points[m_lines[m_faces[fid][0]][0]];
    else
      p0 = m_points[m_lines[-m_faces[fid][0]][1]];

    if(m_faces[fid][1] > 0)
      p1 = m_points[m_lines[m_faces[fid][1]][0]];
    else
      p1 = m_points[m_lines[-m_faces[fid][1]][1]];

    if(m_faces[fid][2] > 0)
      p2 = m_points[m_lines[m_faces[fid][2]][0]];
    else
      p2 = m_points[m_lines[-m_faces[fid][2]][1]];

    auto v0 = p1 - p0;
    auto v1 = p2 - p1;
    n = cross(v0, v1);
  }

  void get_point_tangent(const int eid, const Point p, Vector & t)
  {
    auto & p0 = m_points[m_lines[eid][0]];
    auto & p1 = m_points[m_lines[eid][1]];
    t = p1 - p0;
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
  std::map<int, Cube> & get_volumes()
  {
    return m_volumes;
  }
private:
  std::map<int, Point> m_points;
  std::map<int, Line> m_lines;
  std::map<int, Face> m_faces;
  std::map<int, Cube> m_volumes;
  std::map<int, Sphere> m_spheres;
};

}
}


#endif // end of CubeModel_h

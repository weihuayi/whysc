#ifndef CubeModel_h
#define CubeModel_h

#include <map>
#include<math.h>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class CubeModel
{
public:
  typedef typename GK::Point_3 Point;
  typedef std::array<int, 2> Line;
  typedef std::vector<int> Face;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;
public:
  CubeModel(){ m_num = 0;}

  void add_point(Point & p)
  {
    m_points[m_num] = p;
    m_num++;
  }

  void add_line(int i, int j)
  {
    m_lines[m_num] = Line{i, j};
    m_num++;
  }

  void add_face(std::vector<int> & f)
  {
    m_faces[m_num] = f;
    m_num++;
  }

  void project_to_face(const int fid, Point & p)
  {
      auto f = m_faces[fid]; 
      auto l0 = m_lines[f[0]];
      auto l1 = m_lines[f[1]]
      auto l = std::sqrt(v.squared_length());
      v /= l;
      p = m_center + m_r*v;
  }
  virtual void project_to_edge(const int eid, Point & p) = 0;
  void get_point_normal(const int fid, const Point p, Vector & n)
  {
      n = p - m_center; 
      auto l = std::sqrt(n.squared_length());
      n /= l;
  }

  virtual void get_point_tangent(const int eid, const Point p, Vector & t) = 0;
private:
  int m_num; //点边面的个数之和
  std::map<int, Point> m_points;
  std::map<int, Line> m_lines;
  std::map<int, Face> m_faces;
};

}
}


#endif // end of CubeModel_h

#ifndef SphereModel_h
#define SphereModel_h

#include <map>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class SphereModel
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;
  typedef typename LevelSetFunction::Sphere_3<GK> Sphere;
public:
  Sphere_3():m_center(F(0.0), F(0.0), F(0.0)), m_r(F(1)) {}
  Sphere_3(const F x, const F y, const F z, const F r):m_center(x, y, z), m_r(r){}
  void project_to_face(const int fid, Point & p)
  {
      auto v = p - m_center; 
      auto l = std::sqrt(v.squared_length());
      v /= l;
      p = m_center + m_r*v;
  }
  virtual void project_to_edge(const int eid, Point & p) = 0;
  void get_point_normal(const int fid, const Point p, Vector & n);
  {
      n = p - m_center; 
      auto l = std::sqrt(n.squared_length());
      n /= l;
  }

  virtual void get_point_tangent(const int eid, const Point p, Vector & t) = 0;
private:
    Point m_center;
    F m_r
};

}
}


#endif // end of SphereModel_h
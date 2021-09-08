#ifndef CubeWithSpheresModelHexMesh_h
#define CubeWithSpheresModelHexMesh_h

#include <map>
#include <vector>
#include <math.h>
#include <iostream>

namespace WHYSC {
namespace GeometryModel {

template<class GK>
class CubeWithSpheresModelHexMesh
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;

  typedef typename std::vector<int> Line;
  typedef typename std::vector<int> Face;
  typedef typename std::vector<int> Volume;
  typedef typename GK::Sphere_3 Sphere;

public:
  CubeWithSpheresModelHexMesh(){}

  void project_to_face(const int fid, Point & p)
  {
    double r = std::sqrt(3);
    auto v = p - Point(0, 0, 0);
    v = (r/std::sqrt(v.squared_length()))*v;
    p = Point(0, 0, 0)+v;
  }

  void project_to_edge(const int eid, Point & p){}

  void project_vector_to_face(const int fid, const Point p, Vector & v)
  {
    auto center = Point(0, 0, 0);
    auto n = p - center;
    v = v - dot(v, n)*n/n.squared_length();
  }

  void project_vector_to_edge(const int eid, const Point p, Vector & v) {}
};

}
}


#endif // end of CubeAndSphereModelHexMesh_h

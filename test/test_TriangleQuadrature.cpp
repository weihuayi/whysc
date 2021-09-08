#include <array>
#include <vector>
#include <algorithm>
#include <math.h>

#include "TestMacro.h"
#include "geometry/Point_2.h"
#include "geometry/Vector_2.h"
#include "quadrature/TriangleQuadrature.h"

typedef WHYSC::GeometryObject::Point_2<double> Point;
typedef WHYSC::Quadrature::TriangleQuadrature TriangleQuadrature;

double f(const Point & p)
{
  return p[0];
}

struct Triangle
{
  Point p0, p1, p2;
  double area()
  {
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    return cross(v1, v2)/2;
  }
};

void test_integral(int p = 3)
{
  const double h = 0.5;
  Triangle tri;
  tri.p0 = Point({0, 0});
  tri.p1 = Point({h, 0});
  tri.p2 = Point({0, h});
  
  TriangleQuadrature triQ(p);
  int NQP = triQ.number_of_quadrature_points();

  double val = 0.0;
  for(int i = 0; i < NQP; i++)
  {
    auto & w = triQ.weight[i];
    auto & qpts = triQ.quadrature_point[i];
    auto P = qpts[0]*tri.p0 + qpts[1]*tri.p1 + qpts[2]*tri.p2; 
    val += f(P)*w;
  }
  val *= tri.area();
  std::cout << val << " " << std::pow(h, 3)/6 <<std::endl;
  ASSERT_THROW(val-(std::pow(h, 3)/6)<1e-10);
}

int main(int args, char *argv[])
{
  int q = 4;
  if(args>1)
     q = std::stoi(argv[1]);
  test_integral(q);
}

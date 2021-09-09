
/**
 * \file test_TetrahedronQuadrature.cpp
 * \author Chunyu Chen
 * \date 2021/09/08
 * \brief TetrahedronQuadrature.h 测试文件
 */
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>

#include "TestMacro.h"
#include "geometry/Point_3.h"
#include "geometry/Vector_3.h"
#include "quadrature/TetrahedronQuadrature.h"

typedef WHYSC::GeometryObject::Point_3<double> Point;
typedef WHYSC::Quadrature::TetrahedronQuadrature TetrahedronQuadrature;

/** 
 * \brief 被积函数
 * \param p 函数的参数
 */
double f(const Point & p)
{
  return p[0];//p[0]*p[1];
}

/**
 * \brief 四面体类, 即积分区域
 */
struct Tetrahedron
{
  Point p0, p1, p2, p3; /**< 四面体顶点 */

  /**
   * \brief 四面体面积
   */
  double area()
  {
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    auto v3 = p3 - p0;
    return std::abs(dot(cross(v1, v2), v3)/6);
  }
};

/**
 * \brief 积分测试函数
 */
void test_integral(int p = 3, double h = 0.5)
{
  Tetrahedron tet;
  tet.p0 = Point({0, 0, 0});
  tet.p1 = Point({h, 0, 0});
  tet.p2 = Point({0, h, 0});
  tet.p3 = Point({0, 0, h});
  
  TetrahedronQuadrature tetQ(p);
  int NQP = tetQ.number_of_quadrature_points();

  double val = 0.0;
  for(int i = 0; i < NQP; i++)
  {
    auto w = tetQ.quadrature_weight(i);
    auto & qpts = tetQ.quadrature_point(i);
    auto P = qpts[0]*tet.p0 + qpts[1]*tet.p1 + qpts[2]*tet.p2 + qpts[3]*tet.p3; 
    val += f(P)*w;
  }
  val *= tet.area();
  ASSERT_THROW(std::abs(val-std::pow(h, 5)/12.0)<1e-10);/**< 判断积分是否正确 */
}

int main(int args, char *argv[])
{
  int q = 4;
  double h = 0.5;
  if(args>1 && args<3)
     q = std::stoi(argv[1]);
  else if(args>=3)
  {
     q = std::stoi(argv[1]);
     h = std::stoi(argv[2]);
  }
  test_integral(q, h);
}

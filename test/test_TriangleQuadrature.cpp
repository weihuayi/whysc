
/**
 * \file test_TriangleQuadrature.cpp
 * \author Chunyu Chen
 * \date 2021/09/08
 * \brief TriangleQuadrature.h 测试文件
 */
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

/** 
 * \brief 被积函数
 * \param p 函数的参数
 */
double f(const Point & p)
{
  return p[0];
}

/**
 * \brief 三角形类, 即积分区域
 */
struct Triangle
{
  Point p0, p1, p2; /**< 三角形顶点 */

  /**
   * \brief 三角形面积
   */
  double area()
  {
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    return cross(v1, v2)/2;
  }
};

/**
 * \brief 积分测试函数
 */
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
    auto w = triQ.quadrature_weight(i);
    auto & qpts = triQ.quadrature_point(i);
    auto P = qpts[0]*tri.p0 + qpts[1]*tri.p1 + qpts[2]*tri.p2; 
    val += f(P)*w;
  }
  val *= tri.area();
  std::cout << val << " " << std::pow(h, 3)/6 <<std::endl;
  ASSERT_THROW(val-(std::pow(h, 3)/6)<1e-10);/**< 判断积分是否正确 */
}

int main(int args, char *argv[])
{
  int q = 4;
  if(args>1)
     q = std::stoi(argv[1]);
  test_integral(q);
}

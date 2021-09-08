/**
 * \file test_GaussLegendreQuadrature.cpp
 * \author Chunyu Chen
 * \date 2021/09/08
 * \brief GaussLegendreQuadrature.h 测试文件
 */
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>

#include "TestMacro.h"
#include "geometry/Point_2.h"
#include "geometry/Vector_2.h"
#include "quadrature/GaussLegendreQuadrature.h"

typedef WHYSC::GeometryObject::Point_2<double> Point;
typedef WHYSC::Quadrature::GaussLegendreQuadrature GaussLegendreQuadrature;


/** 
 * \brief 被积函数
 * \param p 函数的参数
 */
double f(const double & x)
{
  return x*x;
}

/**
 * \brief 区间类, 即积分区域
 */
struct Interval
{
  Interval(double a0, double b0): a(a0), b(b0) {}

  double a, b; /**< 区间端点 */

  /**
   * \brief 区间长度
   */
  double area()
  {
    return std::abs(b-a);
  }
};

/**
 * \brief 积分测试函数
 */
void test_integral(int p = 3)
{
  const double h = 0.5;
  Interval I(0, h);
  GaussLegendreQuadrature IQ(p);
  int NQP = IQ.number_of_quadrature_points();

  double val = 0.0;
  for(int i = 0; i < NQP; i++)
  {
    auto w = IQ.quadrature_weight(i);
    auto & qpts = IQ.quadrature_point(i);
    auto P = qpts[0]*I.a + qpts[1]*I.b; 
    val += f(P)*w;
  }
  val *= I.area();
  ASSERT_THROW(val-(std::pow(h, 3)/3)<1e-10);/**< 判断积分是否正确 */
}

int main(int args, char *argv[])
{
  int q = 4;
  if(args>1)
     q = std::stoi(argv[1]);
  test_integral(q);
}

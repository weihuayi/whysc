/**
 * \file GaussLegendreQuadrature.h 
 * 区间 \f$[0, 1]\f$ 上的 Gauss Legendre 积分公式
 * \author Huayi Wei<weihuayi@xtu.edu.cn>
 * \date 2021.09.08
 */
#include <array>
#include <vector>

#include "GaussLegendreQuadratureData.h"

namespace WHYSC{
namespace Quadrature{

/**
 * \class GaussLegendreQuadrature GaussLegendreQuadrature.h
 *
 */
class GaussLegendreQuadrature
{
public:
  /**
   * \brief 一维区间上重心坐标形式的参考点 \f$(\lambda_0, \lambda_1\f$
   * 
   * \f[
   *    \lambda_0 + \lambda_1 = 1
   * \f]
   */
  typedef std::array<double, 2> RNode;

public:

  /**
   * \brief 构造函数
   * q 是积分公式的编号，要求大于等于0, 第 q 个积分公式有 q+1 个积分点
   */
  GaussLegendreQuadrature(int q): m_q(q) 
  {
    int NQ = number_of_quadrature_points();
    m_weight.resize(NQ);
    m_qpoint.resize(NQ);
    
    int s = m_q*(m_q+1)/2;
    for(int i = 0; i < NQ; i++)
    {
      m_weight[i] = GaussLegendreQuadratureData[s+i][1]/2.0;
      m_qpoint[i][0] = GaussLegendreQuadratureData[s+i][0]/2 + 0.5;
      m_qpoint[i][1] = 1.0 - m_qpoint[i][0];
    }
  }

  /**
   * \brief 获取当前积分公式中积分点的个数
   */
  int number_of_quadrature_points() 
  {
    return m_q+1;
  } 

  /**
   * \brief 获取当前积分公式中第 i 个参考积分点
   * \attention 注意是重心坐标形式的积分点，有两个分量 
   */
  const RNode & quadrature_point(int i)
  {
    return m_qpoint[i];
  }

  const double & quadrature_weight(int i)
  {
    return m_weight[i];
  }

  const std::vector<RNode> & quadrature_points()
  {
    return m_qpoint;
  }

  const std::vector<double> & quadrature_weights()
  {
    return m_weight;
  }

private:
  int m_q;
  std::vector<double> m_weight;
  std::vector<RNode>  m_qpoint;
};
}
}

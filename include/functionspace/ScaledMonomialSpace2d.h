#ifndef ScaledMonomialSpace2d_h
#define ScaledMonomialSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>

namespace WHYSC {
namespace FunctionSpace{

template<typename Mesh>
class SMDof2d
{
public:
  typedef std::array<int, 2> Dof;

public:
  SMDof2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p) {}

  int number_of_local_dofs()
  {
    return (m_p+1)*(m_p+2)/2;
  }

  int number_of_dofs()
  {
    auto NC = m_mesh->number_of_cells();
    return NC*number_of_local_dofs();
  }

  void cell_to_dof(int idx, std::vector<int> & cell2dof)
  {
    int ldof = number_of_local_dofs();
    cell2dof.resize(ldof);

    int N = idx*ldof;
    for(int i = 0; i < ldof; i++)
    {
      cell2dof[i] = N+i;
    }
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
};// end of SMDof2d

template<typename Mesh>
class ScaledMonomialSpace2d
{
public:
  typedef SMDof2d<Mesh> SMDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

public:
  ScaledMonomialSpace2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p)
  {
    m_dof = std::make_shared<SMDof>(mesh, p);
    mesh->cell_barycenter(m_cellbarycenter);
    mesh->cell_size(m_cellsize);
  }

  void basis(int cidx, const Node & point, std::vector<double> & phi)
  {
    int ldof = m_dof->number_of_local_dofs();
    phi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    auto xbar = (point - cellbar)/h;

    phi[0] = 1; // 0 次基函数
    if(m_p>0)
    {
      phi[1] = xbar[0];
      phi[2] = xbar[1];// 1 次基函数
      int start = 3; // 第 0 个 j 次基函数的编号
      for(int j = 2; j < m_p+1; j++)
      {
        for(int k = start; k < start+j; k++)
        {
          phi[k] = xbar[0]*phi[k-j]; 
        }
        phi[start+j] = xbar[1]*phi[start-1];
        start += j+1;
      }
    }
  }

  void grad_basis(int cidx, const Node & point, std::vector<Vector> & nphi)
  {
    int ldof = m_dof->number_of_local_dofs();
    nphi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    auto xbar = (point - cellbar)/h;

    std::vector<double> phi; //基函数 
    basis(cidx, point, phi);

    nphi[0][0] = 0.0;
    nphi[0][1] = 0.0; // 0 次基函数的偏导
    if(m_p>0)
    {
      int start = 1; // 第 0 个 j 次基函数的编号
      for(int j = 1; j < m_p+1; j++)
      {
        nphi[start][1] = 0.0; // phi[-1]
        for(int k = start; k < start+j; k++)
        {
          int c = k - start;
          nphi[k][0] = (j-c)*phi[k-j];
          nphi[k+1][1] = (c+1)*phi[k-j];
        }
        nphi[start+j][0] = 0.0;
        start += j+1; //更新 start
      }
    }
  }

  void laplace_basis(int cidx, const Node & point, std::vector<double> & lphi)
  {
    int ldof = m_dof->number_of_local_dofs();
    lphi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    auto xbar = (point - cellbar)/h;

    std::vector<double> phi; //基函数 
    basis(cidx, point, phi);

    lphi[0] = 0.0; // 0 次基函数的偏导
    if(m_p>0)
    {
      lphi[1] = 0.0;
      lphi[2] = 0.0;
      phi[1] = xbar[0];
      phi[2] = xbar[1];
      int start = 3; // 第 0 个 j 次基函数的编号
      for(int j = 2; j < m_p+1; j++)
      {
        for(int k = start; k < start+j-1; k++)
        {
          int c = k - start;
          lphi[k] += (j-c)*(j-c-1)*phi[k-j-j+1];
          lphi[k+2] += (c+2)*(c+1)*phi[k-j-j+1];
        }
        start += j+1; //更新 start
      }
    }
  }

private:
  int m_p;
  std::vector<Node> m_cellbarycenter;
  std::vector<F> m_cellsize;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<SMDof> m_dof;
};// end of ScaledMonomialSpace2d

}//end of FunctionSpace
}//end of WHYSC
#endif // end of ScaledMonomialSpace2d_h

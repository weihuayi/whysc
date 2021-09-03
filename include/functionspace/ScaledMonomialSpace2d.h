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
  SMDof2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p)
  {
    build_dof();
  }

  int number_of_local_dofs()
  {
    return (m_p+1)*(m_p+2)/2;
  }

  int number_of_dofs()
  {
    auto NC = m_mesh->number_of_cells();
    return NC*number_of_local_dofs();
  }

  template<typename Container>
  void cell_to_dof(int idx, Container & cell2dof)
  {
    int ldof = number_of_local_dofs();
    int N = idx*ldof;
    for(int i = 0; i < ldof; i++)
    {
      cell2dof[i] = N+i;
    }
  }

  void build_dof()
  {
    int ldof = number_of_local_dofs(); 
    m_multiIndex.resize(ldof);
    int idx = 0;
    for(int i = 0; i < m_p+1; i++)
    {
      for(int j = 0; j < i+1; j++)
      {
        m_multiIndex[idx][1] = j; 
        m_multiIndex[idx][0] = i-j; 
        idx++;
      }
    }
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  std::vector<Dof> m_multiIndex;
};// end of SMDof2d

template<typename Mesh, typename AK>
class ScaledMonomialSpace2d
{
public:
  typedef SMDof2d<Mesh> SMDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;

  typedef typename AK::Matrix Matrix;
  typedef typename AK::CSRMatrix CSRMatrix;

public:
  ScaledMonomialSpace2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p)
  {
    m_dof = std::make_shared<SMDof>(mesh, p);
    mesh->cell_barycenter(m_cellbarycenter);
    mesh->cell_size(m_cellsize);
  }

  void basic(const std::vector<Node> & point, Matrix & phi, const std::vector<int> & index)
  {
    int ldof = m_dof->number_of_local_dofs();
    int NP = index.size();
    phi.reshape(NP, ldof);
    for(int i = 0; i < NP; i++)
    {
      int idx = index[i];
      const auto & cellbar = m_cellbarycenter[idx];
      auto & h = m_cellsize[idx];
      auto xbar = (point[i] - cellbar)/h;

      phi[i][0] = 1; // 0 次基函数
      if(m_p>0)
      {
        phi[i][1] = xbar[0];
        phi[i][2] = xbar[1];// 1 次基函数
        int start = 3; // 第 0 个 j 次基函数的编号
        for(int j = 2; j < m_p+1; j++)
        {
          for(int k = start; k < start+j; k++)
          {
            phi[i][k] = xbar[0]*phi[i][k-j]; 
          }
          phi[i][start+j] = xbar[1]*phi[i][start-1];
          start += j+1;
        }
      }
    }
  }

  void basic(const std::vector<Node> & point, Matrix & phi)
  {
    int ldof = m_dof->number_of_local_dofs();
    int NC = m_mesh->number_of_cells();
    phi.reshape(NC, ldof);
    for(int i = 0; i < NC; i++)
    {
      const auto & cellbar = m_cellbarycenter[i];
      auto & h = m_cellsize[i];
      auto xbar = (point[i] - cellbar)/h;

      phi[i][0] = 1; // 0 次基函数
      if(m_p>0)
      {
        int start = 1; // 第 0 个 j 次基函数的编号
        for(int j = 1; j < m_p+1; j++)
        {
          for(int k = start; k < start+j; k++)
          {
            phi[i][k] = xbar[0]*phi[i][k-j]; //phi[:-2]
          }
          phi[i][start+j] = xbar[1]*phi[i][start-1]; // phi[-1]
          start += j+1; //更新 start
        }
      }
    }
  }

  void grad_basic(const std::vector<Node> & point, Matrix & phi_x, Matrix & phi_y, 
      const std::vector<int> & index)
  {
    int ldof = m_dof->number_of_local_dofs();
    int NP = index.size();
    phi_x.reshape(NP, ldof);
    phi_y.reshape(NP, ldof);
    for(int i = 0; i < NP; i++)
    {
      int idx = index[i];
      const auto & cellbar = m_cellbarycenter[idx];
      auto & h = m_cellsize[idx];
      auto xbar = (point[i] - cellbar)/h;
      double phi[ldof]; // j-1 次基函数

      phi[0] = 1.0;
      phi_x[i][0] = 0.0;
      phi_y[i][0] = 0.0; // 0 次基函数的偏导
      if(m_p>0)
      {
        int start = 1; // 第 0 个 j 次基函数的编号
        for(int j = 1; j < m_p+1; j++)
        {
          phi_y[i][start] = 0.0; // phi[-1]
          for(int k = start; k < start+j; k++)
          {
            int c = k - start;                 //phi_y[0]
            phi_x[i][k] = (j-c)*phi[k-j];      //phi_x[:-2]
            phi_y[i][k+1] = (c+1)*phi[k-j];    //phi_y[1:]
            phi[k] = xbar[0]*phi[k-j];         //phi[:-2]
          }
          phi_x[i][start+j] = 0.0;               //phi_x[-1] 
          phi[start+j] = xbar[1]*phi[start-1]; //phi[-1]
          start += j+1; //更新 start
        }
      }
    }
  }

  void grad_basic(const std::vector<Node> & point, Matrix & phi_x, Matrix & phi_y)
  {
    int ldof = m_dof->number_of_local_dofs();
    int NC = m_mesh->number_of_cells();
    phi_x.reshape(NC, ldof);
    phi_y.reshape(NC, ldof);
    for(int i = 0; i < NC; i++)
    {
      const auto & cellbar = m_cellbarycenter[i];
      auto & h = m_cellsize[i];
      auto xbar = (point[i] - cellbar)/h;
      double phi[ldof]; // j-1 次基函数

      phi[0] = 1.0;
      phi_x[i][0] = 0.0;
      phi_y[i][0] = 0.0; // 0 次基函数的偏导
      if(m_p>0)
      {
        int start = 1; // 第 0 个 j 次基函数的编号
        for(int j = 1; j < m_p+1; j++)
        {
          phi_y[i][start] = 0; // phi[-1]
          for(int k = start; k < start+j; k++)
          {
            int c = k - start;                 //phi_y[0]
            phi_x[i][k] = (j-c)*phi[k-j];      //phi_x[:-2]
            phi_y[i][k+1] = (c+1)*phi[k-j];    //phi_y[1:]
            phi[k] = xbar[0]*phi[k-j];         //phi[:-2]
          }
          phi_x[i][start+j] = 0;               //phi_x[-1] 
          phi[start+j] = xbar[1]*phi[start-1]; //phi[-1]
          start += j+1; //更新 start
        }
      }
    }
  }

  void laplace_basic(const std::vector<Node> & point, Matrix & lphi)
  {
    int ldof = m_dof->number_of_local_dofs();
    int NC = m_mesh->number_of_cells();
    lphi.reshape(NC, ldof);
    for(int i = 0; i < NC; i++)
    {
      const auto & cellbar = m_cellbarycenter[i];
      auto & h = m_cellsize[i];
      auto xbar = (point[i] - cellbar)/h;
      double phi[ldof]; // j-1 次基函数

      phi[0] = 1.0;
      lphi[i][0] = 0.0; // 0 次基函数的偏导
      if(m_p>0)
      {
        lphi[i][1] = 0.0;
        lphi[i][2] = 0.0;
        phi[1] = xbar[0];
        phi[2] = xbar[1];
        int start = 3; // 第 0 个 j 次基函数的编号
        for(int j = 2; j < m_p+1; j++)
        {
          for(int k = start; k < start+j-1; k++)
          {
            int c = k - start;
            lphi[i][k] += (j-c)*(j-c-1)*phi[k-j-j+1];
            lphi[i][k+2] += (c+2)*(c+1)*phi[k-j-j+1];
            phi[k] = xbar[0]*phi[k-j];
          }
          phi[start+j-1] = xbar[1]*phi[start-2];
          phi[start+j] = xbar[1]*phi[start-1];
          start += j+1; //更新 start
        }
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

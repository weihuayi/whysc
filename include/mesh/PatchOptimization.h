#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <numeric>

namespace WHYSC {
namespace Mesh {


template<typename Mesh, typename ObjectionFunction, typename Model>
class PatchOptimization
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;
  typedef typename Mesh::Toplogy Toplogy;

public:
  PatchOptimization(std::shared_ptr<Mesh> mesh, std::shared_ptr<Model> model):
    m_mesh(mesh), m_model(model), m_objfun(std::make_shared<ObjectionFunction>(mesh))
  {}

  void optimization(int i)
  {
    Vector move;
    
    computerMoveVector(i, move);

    auto & node = m_mesh->node(i);
    node = node + move;
    project(i, node);
  }

  void project(int i, Node & node)
  {
    auto tag = m_mesh->get_node_int_data()["gtag"][i];
    auto dof = m_mesh->get_node_int_data()["gdof"][i];
    if(dof==2)
    {
      m_model->project_to_face(tag, node);
    }
    else if(dof==1)
    {
      m_model->project_to_edge(tag, node);
    }
  }

  template<typename Patch>
  void bar_optimization(int i, Vector & move, Patch & patch)//点周围单元的重心作为点位置
  {
    int GD = m_mesh->geo_dimension();
    for(int j = 0; j < GD; j++)
      move[j] = -m_mesh->node(i)[j];

    int NP = patch.number_of_adj_entities(i);
    for(int k = 0; k < NP; k++)
    {
      Node tmpNode;
      m_mesh->cell_barycenter(patch.adj_entity(k), tmpNode);
      for(int j = 0; j < GD; j++)
      {
        move[j] += tmpNode[j]/NP;
      }
    }
  }

  void computerMoveVector(int i, Vector & move)//给每个点一个方向, 求这个方向上质量的最小值
  {
    Vector v;
    m_objfun->get_grad_value(i, v);
    auto tmpNode = m_mesh->node(i);

    double a = 0;
    double b = 1;
    double c = a + (1 - 0.618)*(b-a);
    double d = a + 0.618*(b-a);

    m_mesh->node(i) = tmpNode + c*v;
    double Qc = m_objfun->get_value(i);
    m_mesh->node(i) = tmpNode + d*v;
    double Qd = m_objfun->get_value(i);

    while(abs(Qc-Qd)>0.0001)// 0.618法
    {
      if(Qc > Qd)
      {
        a = c;
        c = d;
        d = a + 0.618*(b-a);
        Qc = Qd;
        m_mesh->node(i) = tmpNode + d*v;
        Qd = m_objfun->get_value(i);
      }
      else
      {
        b = d;
        d = c;
        c = a + (1 - 0.618)*(b-a);
        Qd = Qc;
        m_mesh->node(i) = tmpNode + c*v;
        Qc = m_objfun->get_value(i);
      }
    }
    move = c*v;
    m_mesh->node(i) = tmpNode;
  }

  template<typename Patch>
  double min_len(Patch & patch)
  {
    auto & cell = m_mesh->cells();
    auto & node = m_mesh->nodes();
    int NP = patch.number_of_adj_entities();
    double L = 100000.0;

    for(int i = 0; i < NP; i++)
    {
      auto & cell = m_mesh->cell(patch.adj_entity(i));
      auto j = patch.adj_local_index(i);
      auto current = cell[j];

      auto v1 = node[cell[(j+1)%4]] - node[current];
      double L1 = std::sqrt(v1.squared_length());
      if(L>L1)
        L=L1;
    }
    return L;
  }

  std::shared_ptr<Mesh> get_mesh()
  {
    return m_mesh;
  }
  
private:
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<Model> m_model;
  std::shared_ptr<ObjectionFunction> m_objfun;
};

} // end of namespace Mesh

} // end of namespace WHYSC

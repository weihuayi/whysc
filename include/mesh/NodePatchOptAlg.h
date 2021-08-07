#ifndef NodePatchOptAlg_h
#define NodePatchOptAlg_h

#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <numeric>

namespace WHYSC {
namespace Mesh {


template<typename Mesh, typename ObjectFunction, typename Model>
class NodePatchOptAlg
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;
  typedef typename Mesh::Toplogy Toplogy;

public:
  NodePatchOptAlg(std::shared_ptr<Mesh> mesh, std::shared_ptr<Model> model):
    m_mesh(mesh), m_model(model), m_objfun(std::make_shared<ObjectFunction>(mesh))
  {}

  void optimization(int i)
  {
    auto move = computerMoveVector(i);

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

  Vector computerMoveVector(int i)//给每个点一个方向, 求这个方向上质量的最小值
  {
    m_objfun->set_patch(i); //设置当前 patch
    auto v = m_objfun->gradient();

    double a = 0;
    double b = 1;
    double c = a + (1 - 0.618)*(b-a);
    double d = a + 0.618*(b-a);

    double Qc = m_objfun->value(m_mesh->node(i)+c*v);
    double Qd = m_objfun->value(m_mesh->node(i)+d*v);

    while(abs(Qc-Qd)>0.0001)// 0.618法
    {
      if(Qc > Qd)
      {
        a = c;
        c = d;
        d = a + 0.618*(b-a);
        Qc = Qd;
        Qd = m_objfun->value(m_mesh->node(i)+d*v);
      }
      else
      {
        b = d;
        d = c;
        c = a + (1 - 0.618)*(b-a);
        Qd = Qc;
        Qc = m_objfun->value(m_mesh->node(i)+c*v);
      }
    }
    if(m_mesh->id()==1)
      std::cout<<m_objfun->value(m_mesh->node(i)) <<  " " << m_objfun->value(m_mesh->node(i)+c*v) <<std::endl;
    return c*v;
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
  std::shared_ptr<ObjectFunction> m_objfun;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of NodePatchOptAlg_h

#ifndef NodePatchOptAlg_h
#define NodePatchOptAlg_h

#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <numeric>
#include <functional>

namespace WHYSC {
namespace Mesh {


template<typename Mesh, typename Model>
class NodePatchOptAlg
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;
  typedef typename Mesh::Toplogy Toplogy;

public:
  NodePatchOptAlg(std::shared_ptr<Mesh> mesh, std::shared_ptr<Model> model):
    m_mesh(mesh), m_model(model){}

  template<typename ObjectFunction>
  double optimization(ObjectFunction & objfun)
  {
    int i = objfun.get_patch()->id();
    auto & node = m_mesh->node(i);
    auto v = objfun.direction();
    preprocess(i, node, v);

    auto alpha = computer_step(objfun, node, v);

    node = node + alpha*v;
    proprocess(i, node);
    return alpha*std::sqrt(v.squared_length()); 
  }

  void preprocess(int i, const Node & node, Vector & v)
  {
    auto tag = m_mesh->get_node_int_data()["gtag"][i];
    auto dof = m_mesh->get_node_int_data()["gdof"][i];
    if(dof==2)
    {
      m_model->project_vector_to_face(tag, node, v);
    }
    else if(dof==1)
    {
      m_model->project_vector_to_edge(tag, node, v);
    }
  }

  void proprocess(int i, Node & node)
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

  template<typename ObjectFunction>
  double computer_step(ObjectFunction objfun, Node & node, Vector & v)
  {
    double a = 0;
    double b = 1;
    double c = a + (1 - 0.618)*(b-a);
    double d = a + 0.618*(b-a);

    double Qc = objfun.value(node+c*v);
    double Qd = objfun.value(node+d*v);
    while(abs(Qc-Qd)>0.0001 || abs(a-b) > 1e-8)// 0.618法
    {
      if(Qc > Qd)
      {
        a = c;
        c = d;
        d = a + 0.618*(b-a);
        Qc = Qd;
        Qd = objfun.value(node+d*v);
      }
      else
      {
        b = d;
        d = c;
        c = a + (1 - 0.618)*(b-a);
        Qd = Qc;
        Qc = objfun.value(node+c*v);
      }

      if(abs(a-b)<1e-16)
        break;
    }
    return c;
  }

private:
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<Model> m_model;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of NodePatchOptAlg_h

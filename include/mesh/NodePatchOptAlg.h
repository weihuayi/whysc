#ifndef NodePatchOptAlg_h
#define NodePatchOptAlg_h

#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <numeric>
#include <functional>

#include "OptimizationAlg.h"

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
    m_mesh(mesh), m_model(model){}

  double optimization(ObjectFunction & objfun)
  {
    int i = objfun.get_patch()->id();
    auto & node = m_mesh->node(i);
    auto v = -1.0*objfun.gradient();
    preprocess(i, node, v);

    std::function<double(double)> F = [&objfun, v](double k){return objfun.value(k*v);}; 
    auto k = OptimizationAlg::line_search(F);

    node = node + k*v;
    proprocess(i, node);
    return k*std::sqrt(v.squared_length()); 
  }

  void preprocess(int i, const Node & node, Vector & v)
  {
    auto l = std::sqrt(v.squared_length());
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
    v = l*v/std::sqrt(v.squared_length());
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

private:
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<Model> m_model;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of NodePatchOptAlg_h

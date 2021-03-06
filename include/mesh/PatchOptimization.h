#include <memory>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <numeric>

namespace WHYSC {
namespace Mesh {


template<typename Mesh, typename MeshQuality, typename Model>
class PatchOptimization
{
public:
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;
  typedef typename Mesh::Toplogy Toplogy;

public:
  PatchOptimization(std::shared_ptr<Mesh> mesh, std::shared_ptr<Model> model)
  {
    m_mesh = mesh;
    m_model = model;
    m_quality = std::make_shared<MeshQuality>(mesh);
    //mesh->node_to_cell(m_top);
    //mesh->node_to_node(m_n2n);
  }

  template<typename Patch>
  void optimization(int i, Node & node, Patch & patch)
  {
    Vector move;
    computerMoveVector(i, move, patch);

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

  template<typename Patch>
  void computerMoveVector(int i, Vector & move, Patch & patch)//给每个点一个方向, 求这个方向上质量的最小值
  {
    int NP = patch.number_of_adj_entities();
    std::vector<int> cidx(NP);
    Vector v = {0};
    double w;
    for(int k = 0; k < NP; k++)
    {
      cidx[k] = patch.adj_entity(k);
      Vector tmpVector;
      double tmpW;
      m_quality->nabla_quality(patch.adj_entity(k), patch.adj_local_index(k), tmpVector, tmpW);
      v = v + tmpVector;
      w = w + tmpW;
    }
    v = v/w;
    //v = v/std::sqrt(v.squared_length());
    auto tmpNode = m_mesh->node(i);

    double a = 0;
    double b = 1; //min_edge_length(i);
    double c = a + (1 - 0.618)*(b-a);
    double d = a + 0.618*(b-a);

    m_mesh->node(i) = tmpNode + c*v;
    double Qc = m_quality->patch_quality(cidx);
    m_mesh->node(i) = tmpNode + d*v;
    double Qd = m_quality->patch_quality(cidx);

    while(abs(Qc-Qd)>0.0001)// 0.618法
    {
      if(Qc > Qd)
      {
        a = c;
        c = d;
        d = a + 0.618*(b-a);
        Qc = Qd;
        m_mesh->node(i) = tmpNode + d*v;
        Qd = m_quality->patch_quality(cidx);
      }
      else
      {
        b = d;
        d = c;
        c = a + (1 - 0.618)*(b-a);
        Qd = Qc;
        m_mesh->node(i) = tmpNode + c*v;
        Qc = m_quality->patch_quality(cidx);
      }
    }
    move = c*v;
    m_mesh->node(i) = tmpNode;
  }

/*
  double min_edge_length(int k)//点 k 周围的边的长度最小值
  {
    auto & loc = m_n2n.locations();
    auto & nei = m_n2n.neighbors();
    double l = 10000000;
    for(int i = loc[k]; i < loc[k+1]; i++)
    {
      auto v = m_mesh->node(nei[i]) - m_mesh->node(k);
      auto l0 = std::sqrt(v.squared_length());
      if(l > l0)
        l = l0;
    }
    return l;
  }

  int local_index(int i, int j)
  {
    auto & loc = m_top.locations();
    auto & locidx = m_top.local_indices();
    return locidx[loc[i]+j];
  }

  int patch_to_cell(int i, int j)
  {
    auto & loc = m_top.locations();
    auto & nei = m_top.neighbors();
    return nei[loc[i]+j];
  }

  int patch_size(int i)
  {
    auto & loc = m_top.locations();
    return loc[i+1] -loc[i];
  }

*/
  std::shared_ptr<Mesh> get_mesh()
  {
    return m_mesh;
  }
  
private:
  //Toplogy m_top;
  //Toplogy m_n2n;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<MeshQuality> m_quality;
  std::shared_ptr<Model> m_model;
};

} // end of namespace Mesh

} // end of namespace WHYSC

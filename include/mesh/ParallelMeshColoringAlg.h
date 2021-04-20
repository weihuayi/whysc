#include <memory>
#include <vector>
#include <mpi.h>


namespace WHYSC {
namespace Mesh {

template<typename PMesh, typename Communication>
class ParallelMeshColoringAlg
{
public:
  typedef typename PMesh::Toplogy Toplogy;

public:
  ParallelMeshColoringAlg(std::shared_ptr<PMesh> mesh, MPI_Comm comm)
  {
    m_mesh = mesh;
    m_comm = comm;
  }

  void mesh_coloring(std::vector<int> & color)
  {
    auto mesh = get_mesh();
    auto NN = mesh->number_of_nodes();
    auto NE = mesh->number_of_edges();
    auto & pds = mesh->parallel_data_structure();

    Communication communication(mesh, m_comm);
    auto & isGhostNode = communication.get_ghost_node();

    std::vector<int> randVal(NN);
    std::vector<bool> isMin(NN, true);

    std::list<int> edges; //没有被删除的边
    std::list<int> nColored; //没有被染色的点

    for(int i = 0; i < NN; i++)
    {
      if(!isGhostNode[i])
        nColored.push_back(i);
    }

    for(int i=0; i < NE; i++)
      edges.push_back(i);

    int cmax = 0;
    int tnum = 1;
    while(tnum > 0)
    {
      cmax++;
      bool isMin[NN];
      for(auto idx: nColored) // 只对没有染色的点进行循环
      {
        isMin[idx] = true; //初始所有未染色点都是最小点
        randVal[idx] = rand();
      }

      communication.communicate(randVal); // 通信影像节点上的随机值

      for(auto it = edges.begin(); it != edges.end();)
      {
        auto e = mesh->edge(*it);
        //判断*i是否为局部最小data值
        if(color[e[0]]!=0 || color[e[1]]!=0)
        {
          it = edges.erase(it); // 删除当前的边编号, 并移动到下一个边
        }

        else
        {
          if(randVal[e[0]] < randVal[e[1]])
          {
            isMin[e[1]] = false;
          }
          else
          {
            isMin[e[0]] = false;
          }
          it++;
        }
      }

      for(auto it = nColored.begin(); it != nColored.end();)
      {
        if(isMin[*it])
        {
          color[*it] = cmax;
          it = nColored.erase(it);
        }
        else
        {
          it++;
        }
      }

      communication.communicate(color); // 通信影像节点上的颜色值

      int lnum = nColored.size();
      MPI_Allreduce(&lnum, &tnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  void color_test(std::vector<int> & color)
  {
    //检验染色是否成功
    int a=0;
    auto mesh = get_mesh();
    int NC = mesh->number_of_cells();
    int NE = mesh->number_of_edges();
    int NN = mesh->number_of_nodes();

    Toplogy node2node;
    mesh->node_to_node(node2node);
    auto & nei = node2node.neighbors();
    auto & loc = node2node.locations();

    for(int i = 0; i < NN; i++)
    {
      for(int j = loc[i]; j < loc[i+1]; j++)
      {
        if(color[i]==color[nei[j]])
        {
          a++;
        }
      }
    }
    if(a==0)
    {
      std::cout<< "染色成功  单元数 " << NC << " 节点数 " << NN << " 边数 " << NE <<endl; 
    }
    else
    {
      std::cout<< "染色失败  单元数 " << NC << " 节点数 " << NN << " 边数 " << NE <<endl; 
    }
  }


  std::shared_ptr<PMesh> get_mesh()
  {
    return m_mesh;
  }

private:
  MPI_Comm m_comm;
  std::shared_ptr<PMesh> m_mesh;
};


} // end of namespace Mesh

} // end of namespace WHYSC

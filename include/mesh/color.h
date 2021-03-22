#include<stdio.h>
#include<iostream>
#include<list>
#include<vector>
#include<set>


template<typename Mesh>
void mesh_coloring(Mesh & mesh, int * color)
{
  typedef typename Mesh::Toplogy Toplogy;

  auto NN = mesh.number_of_nodes();
  auto NE = mesh.number_of_edges();

  int c = 0;
  std::vector<int> randVal(NN);
  std::vector<bool> isMin(NN, true);

  std::list<int> edges; //没有被删除的边
  std::list<int> nColored; //没有被染色的点

  for(int i=0; i < NE; i++)
  {
    edges.push_back(i);
  }

  for(int i=0; i < NN; i++)
  {
    nColored.push_back(i);
  }

  //染色
  while(nColored.size()>0)
  {
    c++;
    // 只对没有染色的点进行循环
    for(const int& index: nColored)
    {
      isMin[index] = true; //初始所有未染色点都是最小点
      randVal[index] = rand();
    }

    for(auto it = edges.begin(); it != edges.end();)
    {
      auto e = mesh.edge(*it);
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
        color[*it] = c;
        it = nColored.erase(it);
      }
      else
      {
        it++;
      }
    }

  }

  //降低颜色数量
  int cMax = c;
  Toplogy node2node;
  mesh.node_to_node(node2node);

  // 尝试用新的 API
  auto & nei = node2node.neighbors();
  auto & loc = node2node.locations();

  while(true)
  {
    for(int i = 0; i < NN; i++)
    {
      std::set<int> com;
      for(int k=1; k <= cMax; k++)
      {
        com.insert(k);
      }

      for(int j = loc[i]; j < loc[i+1]; j++)
      {
          com.erase(color[nei[j]]);
      }
      color[i] = *std::min_element(com.begin(), com.end());
    }

    bool flag = true;
    for(int i = 0; i < NN; i++)
    {
        if( color[i] == cMax )
        {
          flag == false;
          break;
        }
    }

    if(flag)
    {
      cMax -= 1;
    }
    else
    {
      break;
    }
  }
};


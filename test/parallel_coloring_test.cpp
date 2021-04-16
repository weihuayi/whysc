#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PMesh;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::Mesh::EntityOverlap<int> EntityOverlap;

template<typename I>
void communication(std::shared_ptr<PMesh> mesh, std::vector<I> & data, int dim)
{
  // 通信影像节点上的随机值
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto LNN = mesh->number_of_local_nodes();
  auto & pds = mesh->parallel_data_structure();
  for(auto map : pds)
  {
    auto & target = map.first; 
    auto & meshOverlap = map.second; 
    auto & overlap = meshOverlap.entity_overlap(dim);

    auto & locid = overlap.loc_index();
    auto & adjid = overlap.adj_index();

    int N = adjid.size();
    int adjData[2*N];

    for(int j = 0; j < N; j++)
    {
      adjData[j*2] = -1;
      if(locid[j]<LNN)//只发送自己的数据
      {
        adjData[j*2] = adjid[j];
        adjData[j*2+1] = data[locid[j]];
        if(rank == 1 & target == 0)
          std::cout<< "1->0 " << adjData[2*j] << " " << adjData[2*j+1] <<std::endl;
      }
    }
    MPI_Send(adjData, N*2, MPI_INT, target, 1, MPI_COMM_WORLD);
  }//发送数据完成

  for(auto map : pds)
  {
    auto & target = map.first; 
    auto & meshOverlap = map.second; 
    auto & overlap = meshOverlap.entity_overlap(dim);
    auto & adjid = overlap.adj_index();

    int N = adjid.size();
    int locData[2*N];

    MPI_Recv(locData, N*2, MPI_INT, target, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int k = 0; k < N; k++)
    {
      if(locData[2*k]>=LNN)//只接收别人的数据
      {
        data[locData[2*k]] = locData[k*2+1];//填充影像节点数据
        if(rank == 0 & target == 1)
          std::cout<< "0<-1 " << locData[2*k] << " " << locData[2*k+1] <<std::endl;
      }
    }
  }//接收数据完成
}


void mesh_coloring(std::shared_ptr<PMesh> mesh, std::vector<int> & color)
{
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto NN = mesh->number_of_nodes();
  auto NE = mesh->number_of_edges();
  auto LNN = mesh->number_of_local_nodes();
  auto & pds = mesh->parallel_data_structure();

  std::vector<int> randVal(NN);
  std::vector<bool> isMin(NN, true);

  std::list<int> edges; //没有被删除的边
  std::list<int> nColored; //没有被染色的点

  for(int i = 0; i < LNN; i++)
    nColored.push_back(i);

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

    //if(rank==1)
    //{
    //  std::cout<< "hahahah" << randVal[21] <<std::endl;
    //}
    //communication(mesh, randVal, 0); // 通信影像节点上的随机值
    //if(rank==1)
    //{
    //  std::cout<< "hahahah1" << randVal[21] <<std::endl;
    //}

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
      //if(rank==1)
      //{
      //  std::cout<< "爷还在" << e[0] << " " << isMin[e[0]] << " " << isMin[e[1]]<<std::endl;
      //  std::cout<< randVal[21] <<std::endl;
      //}
      //if(rank==2)
      //{
      //  std::cout<< "y还在" << randVal[20];
      //}
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

    communication(mesh, color, 0); // 通信影像节点上的随机值
    if(rank==1)
    {
      //for(auto i : nColored)
        //std::cout<< i <<std::endl;
        //std::cout<< LNN <<std::endl;
    }

    int lnum = nColored.size();
    MPI_Allreduce(&lnum, &tnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //std::cout<< "lnum = " << lnum << " " << rank <<std::endl;
  }
}


int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.build_mesh();

  auto NN = mesh->number_of_nodes();
  std::vector<int> color(NN);

  mesh_coloring(mesh, color);

  //检验染色是否成功
  int a=0;
  int NC = mesh->number_of_cells();
  int NE = mesh->number_of_edges();
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

  MPI_Finalize();
  return 0;
}

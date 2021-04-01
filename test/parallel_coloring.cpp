#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/VTKMeshReader.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PMesh;
typedef PMesh::Cell Cell;
typedef PMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;

template<typename I>
void communication(PMesh & mesh, std::vector<I> & data, int tag)
{
  // 通信影像节点上的随机值
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto & pds = mesh.parallel_data_structure();
  auto & gid = mesh.node_global_id();
  auto & ng2l = mesh.node_global_to_local_id();
  auto & num = mesh.number_of_nodes_in_process();

  //MPI_Barrier(MPI_COMM_WORLD);
  for(auto it = pds.begin(); it != pds.end(); it++)
  {
    auto key = it->first; 
    auto value = it->second; 
    if(key != rank)
    {
      int N = value.size();
      int gid_data[2*N];
      int j = 0;
      for(auto v : value)
      {
        gid_data[j*2] = gid[v];
        gid_data[j*2+1] = data[v];
        j++;
      }

      MPI_Send(&N, 1, MPI_INT, key, tag*10, MPI_COMM_WORLD);
      MPI_Send(gid_data, N*2, MPI_INT, key, tag*10+1, MPI_COMM_WORLD);
      //std::cout<< "myrank = " << rank << " 发给 " << key << " 完成 "<<std::endl; 
    }
  }//发送数据完成
  //MPI_Barrier(MPI_COMM_WORLD);

  for(int j = 0; j < nprocs; j++)
  {
    if(j != rank)
    {
      int N;
      MPI_Recv(&N, 1, MPI_INT, j, tag*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int gid_data[N*2];
      MPI_Recv(gid_data, N*2, MPI_INT, j, tag*10+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //std::cout<< "myrank = " << rank << " 接收 " << j << " 完成 "<<std::endl; 
      for(int k = 0; k < N; k++)
      {
        data[ng2l[gid_data[k*2]]] = gid_data[k*2+1];//填充影像节点数据
      }
    }
  }//接收数据完成
}

void mesh_coloring(PMesh & mesh, std::vector<int> & color)
{
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto NN = mesh.number_of_nodes();
  auto NE = mesh.number_of_edges();

  color.resize(NN);

  std::vector<int> randVal(NN);
  std::vector<bool> isMin(NN, true);

  std::list<int> edges; //没有被删除的边
  std::list<int> nColored; //没有被染色的点
  auto & pds = mesh.parallel_data_structure();

  std::cout<< mesh.number_of_local_nodes() << std::endl;
  for(int i = 0; i < mesh.number_of_local_nodes(); i++)
    nColored.push_back(i);

  for(int i=0; i < NE; i++)
    edges.push_back(i);

  int cmax = 0;
  int tnum = 1;
  while(tnum > 0)
  {
    cmax++;
    bool isMin[NN];
     // 只对没有染色的点进行循环
    for(const int& idx: nColored)
    {
      isMin[idx] = true; //初始所有未染色点都是最小点
      randVal[idx] = rand();
    }

    // 通信影像节点上的随机值
    std::cout<< "randVal 开始通信"<<std::endl;
    communication(mesh, randVal, 0);
    std::cout<< "randVal 通信完成"<<std::endl;

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
        color[*it] = cmax;
        it = nColored.erase(it);
      }
      else
      {
        it++;
      }
    }
    // 通信影像节点上的随机值
    communication(mesh, color, 1);

    int lnum = nColored.size();
    MPI_Allreduce(&lnum, &tnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::cout<< "tnum = " << tnum <<std::endl;
  }
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream ss;
  ss << argv[1] << "_" << rank << ".vtu";
  std::string fname = ss.str();

  std::cout << "process " << rank << " of " << nprocs << " process." 
    << "The file name is " << fname <<  std::endl;

  PMesh mesh(rank);
  Reader reader(&mesh);
  reader.read(ss.str());

  auto & npid = mesh.node_process_id();
  reader.get_node_data("nid", npid);

  auto & gid = mesh.node_global_id();
  reader.get_node_data("gid", gid);

  mesh.construct_parallel_data_structure();
  std::vector<int> color; 
  std::cout<< "here11" <<std::endl;
  mesh_coloring(mesh, color);

  MPI_Finalize();
  return 0;
}

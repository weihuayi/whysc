#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/C6H6.h"
#include "geometry/C6.h"
#include "geometry/Geometry_kernel.h"
#include "geometry/CubeWithSpheresModelHexMesh.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshOptAlg.h"
#include "mesh/HexJacobiQuality.h"
#include "mesh/HexPositiveJacobiQuality.h"
#include "mesh/SumNodePatchObjectFunction.h"
#include "mesh/MaxNodePatchObjectFunction.h"
#include "mesh/MixNodePatchObjectFunction.h"
#include "mesh/MeshFactory.h"
#include "Python.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
//typedef WHYSC::GeometryModel::CubeWithSpheresModelHexMesh<GK> Model;
//typedef WHYSC::GeometryModel::C6H6<GK> Model;
typedef WHYSC::GeometryModel::C6<GK> Model;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, HexMesh> PMesh;
//typedef WHYSC::Mesh::HexJacobiQuality<PMesh> CellQuality;
typedef WHYSC::Mesh::HexPositiveJacobiQuality<PMesh> CellQuality;
typedef WHYSC::Mesh::SumNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
//typedef WHYSC::Mesh::MaxNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
//typedef WHYSC::Mesh::MixNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef WHYSC::Mesh::ParallelMeshOptAlg<PMesh, ObjectFunction, Model> PMeshOpt;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;

template<typename I>
void plot(std::vector<I> & data0, std::vector<I> & data1)
{
  int N0 = data0.size();
  int N1 = data1.size();

  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('../example')");

  PyObject* pModule = PyImport_ImportModule("plot");
  PyObject* pFunc = PyObject_GetAttrString(pModule, "Histogram_plot");//要运行的函数

  PyObject* plist0 = PyList_New(N0);//函数的参数是一个list
  PyObject* ptuple0 = PyTuple_New(1);//把参数用 tuple 装起来
  PyObject* plist1 = PyList_New(N1);//函数的参数是一个list
  PyObject* ptuple1 = PyTuple_New(1);//把参数用 tuple 装起来

  for(int i = 0; i < N0; i++)
  {
    PyObject*  pra = Py_BuildValue("d", data0[i]);
    PyList_SetItem(plist0, i, pra);
  }

  for(int i = 0; i < N1; i++)
  {
    PyObject*  pra = Py_BuildValue("d", data1[i]);
    PyList_SetItem(plist1, i, pra);
  }

  PyTuple_SetItem(ptuple0, 0, plist0);
  PyTuple_SetItem(ptuple1, 0, plist1);

	PyObject_CallObject(pFunc, ptuple0);//运行函数
	PyObject_CallObject(pFunc, ptuple1);//运行函数
	Py_Finalize();
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.get_mesh();

  auto quad = std::make_shared<Model>();

  auto NC = mesh->number_of_cells();
  std::vector<double> cellQualityInit(NC);
  std::vector<double> cellQualityOpt(NC);

  CellQuality mq(mesh);
  mq.quality_of_mesh(cellQualityInit);

  PMeshOpt optAlg(mesh, quad, MPI_COMM_WORLD);
  clock_t S = clock();
  optAlg.optimization(1e-4, 100);//优化
  clock_t E = clock();
  double runtime = double (E-S)/CLOCKS_PER_SEC;

  mq.quality_of_mesh(cellQualityOpt);

  int NAC = 0;
  MPI_Allreduce(&NC, &NAC, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(mesh->id() == 0)
  {
    //plot(cellQualityInit, cellQualityOpt);
    double i_min_q = 100000000.0;
    double i_max_q = 0.0;
    double i_ave_q = 0.0;
    double o_min_q = 100000000.0;
    double o_max_q = 0.0;
    double o_ave_q = 0.0;

    for(int i = 0; i < NC; i++)
    {
      i_min_q = i_min_q>cellQualityInit[i]? cellQualityInit[i]:i_min_q;
      i_max_q = i_max_q>cellQualityInit[i]? i_max_q:cellQualityInit[i];
      i_ave_q += cellQualityInit[i]/NC;

      o_min_q = o_min_q>cellQualityOpt[i]? cellQualityOpt[i]:o_min_q;
      o_max_q = o_max_q>cellQualityOpt[i]? o_max_q:cellQualityOpt[i];
      o_ave_q += cellQualityOpt[i]/NC;
    }
    std::cout<< "总运行时间:" << runtime << std::endl;
    std::cout<< "初始网格单元最差质量为: " << i_ave_q <<std::endl; 
    std::cout<< "初始网格单元最优质量为: " << i_min_q <<std::endl; 
    std::cout<< "初始网格单元平均质量为: " << i_max_q <<std::endl; 
    std::cout<< "优化网格单元最差质量为: " << o_ave_q <<std::endl; 
    std::cout<< "优化网格单元最优质量为: " << o_min_q <<std::endl; 
    std::cout<< "优化网格单元平均质量为: " << o_max_q <<std::endl; 

    std::ofstream wfile;
    wfile.open("hex_test.txt", ios::app);
    wfile<< argv[1] << "\n\n";
    wfile<< "总运行时间为: " << runtime << "\n";
    wfile<< "单元数目为: " << NAC << "\n";
    wfile<< "进程数为: " << nprocs << "\n";
    wfile<< "\n";
    wfile<< "初始网格单元最差质量为: " << i_ave_q << "\n";
    wfile<< "初始网格单元最优质量为: " << i_min_q << "\n";
    wfile<< "初始网格单元平均质量为: " << i_max_q << "\n";
    wfile<< "\n";
    wfile<< "优化网格单元最差质量为: " << o_ave_q << "\n";
    wfile<< "优化网格单元最优质量为: " << o_min_q << "\n";
    wfile<< "优化网格单元平均质量为: " << o_max_q << "\n";
    wfile<< "\n\n\n\n";
    wfile.close();
  }

  std::stringstream ss;
  ss << "opt_"<< mesh->id() << ".vtu";


  std::vector<int> qInit(NC);
  std::vector<int> qOpt(NC);
  for(int i = 0; i < NC; i++)
  {
    qInit[i] = (int)(cellQualityInit[i]*1000);
    qOpt[i] = (int)(cellQualityOpt[i]*1000);
  }

  Writer writer;
  writer.set_points(*mesh);
  writer.set_cells(*mesh);
  writer.set_point_data(mesh->get_node_int_data()["nid"], 1, "nid");
  writer.set_cell_data(qOpt, 1, "q_opt");
  writer.set_cell_data(qInit, 1, "q_init");
  writer.write(ss.str());

  MPI_Finalize();
  return 0;
}

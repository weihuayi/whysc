#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "geometry/CubeWithSpheresModel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshColoringAlg.h"
#include "mesh/ParallelMeshOptimization.h"
#include "mesh/TetRadiusRatioQualityFunction.h"
#include "mesh/MeshFactory.h"
#include "Python.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::GeometryModel::CubeWithSpheresModel<GK> Model;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TetMesh> PMesh;
typedef WHYSC::Mesh::TetRadiusRatioQualityFunction<PMesh> ObjectionFunction;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef WHYSC::Mesh::ParallelMeshColoringAlg<PMesh> PCA;
typedef WHYSC::Mesh::ParallelMeshOptimization<PMesh, ObjectionFunction, Model> PMeshOpt;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
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

  auto cube = std::make_shared<Model>(2);

  auto NC = mesh->number_of_cells();
  std::vector<double> cellQualityInit(NC);
  std::vector<double> cellQualityOpt(NC);
  ObjectionFunction mq(mesh);
  for(int i = 0; i < NC; i++)
  {
    cellQualityInit[i] = 1/mq.quality_of_cell(i);
  }

  std::cout<< "开始染色..." <<std::endl;
  PCA colorAlg(mesh, MPI_COMM_WORLD);
  colorAlg.coloring();//染色
  colorAlg.color_test();//染色测试

  PMeshOpt optAlg(mesh, cube, MPI_COMM_WORLD);
  for(int i = 0; i < 50; i++)
  {
    std::cout<< "正在优化第 " << i+1 << " 次" <<std::endl;
    optAlg.mesh_optimization();//优化
  }

  for(int i = 0; i < NC; i++)
  {
    cellQualityOpt[i] = 1/mq.quality_of_cell(i);
  }

  if(mesh->id() == 1)
  {
    plot(cellQualityInit, cellQualityOpt);
  }

  std::stringstream ss;
  ss << "opt_"<< mesh->id() << ".vtu";

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();
  writer.write(ss.str());

  MPI_Finalize();
  return 0;
}

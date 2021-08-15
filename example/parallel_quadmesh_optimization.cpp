#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "geometry/RectangleWithHole.h"
#include "mesh/QuadMesh.h"
#include "mesh/ParallelMeshNew.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshOptAlg.h"
#include "mesh/QuadJacobiQuality.h"
#include "mesh/SumNodePatchObjectFunction.h"
#include "mesh/MaxNodePatchObjectFunction.h"
#include "mesh/MeshFactory.h"
#include "Python.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::GeometryModel::RectangleWithHole<GK> Model;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, QMesh> PMesh;
typedef WHYSC::Mesh::QuadJacobiQuality<PMesh> CellQuality;
typedef WHYSC::Mesh::SumNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
//typedef WHYSC::Mesh::MaxNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
typedef WHYSC::Mesh::ParallelMesher<PMesh> PMesher;
typedef WHYSC::Mesh::ParallelMeshOptAlg<PMesh, ObjectFunction, Model> PMeshOpt;
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

  auto quad = std::make_shared<Model>();

  auto NC = mesh->number_of_cells();
  std::vector<double> cellQualityInit(NC);
  std::vector<double> cellQualityOpt(NC);

  CellQuality mq(mesh);
  mq.quality_of_mesh(cellQualityInit);

  PMeshOpt optAlg(mesh, quad, MPI_COMM_WORLD);
  optAlg.optimization(1e-4, 30);//优化

  mq.quality_of_mesh(cellQualityOpt);

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

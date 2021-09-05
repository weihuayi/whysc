#include <sstream> 
#include <cmath> 

#include <fstream>
#include <metis.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/QuadMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef WHYSC::Mesh::QuadMesh<GK, Node, Vector> QuadMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, TriMesh> PTMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, QuadMesh> PQMesh;

typedef WHYSC::Mesh::VTKMeshReader<PTMesh> TReader;
typedef WHYSC::Mesh::VTKMeshWriter<PTMesh> TWriter;
typedef WHYSC::Mesh::VTKMeshReader<PQMesh> QReader;
typedef WHYSC::Mesh::VTKMeshWriter<PQMesh> QWriter;
typedef WHYSC::Mesh::MeshFactory MF;

const double pi = acos(-1.0);

int main() 
{
  /*
  PTMesh mesh(1);
  MF::one_triangle_mesh(mesh);
  mesh.uniform_refine(9);

  std::vector<PTMesh> submeshes;
  */

  PQMesh mesh(1);

  MF::three_quad_mesh(mesh);
  mesh.uniform_refine(4);
  std::vector<PQMesh> submeshes;

  QWriter writer(&mesh);
  writer.set_points();
  writer.set_cells();
  writer.write("init_mesh.vtu");

  int nparts = 4;
  submeshes.resize(nparts);

  auto & cells = mesh.cells();
  auto & nodes = mesh.nodes();
  auto NN = mesh.number_of_nodes();

 //找到网格重心
  Node node_bar;
  for(auto & node : nodes)
  {
    node_bar[0] += node[0]/NN;
    node_bar[1] += node[1]/NN;
    node_bar[2] += node[2]/NN;
  }
  std::vector<int> nid(NN);
  std::vector<int> cid(1);
  for(int k = 0; k < NN; k++)
  {
    auto node = nodes[k];
    double theta = atan2(node[1] - node_bar[1], node[0] - node_bar[0]) + pi;
    for(int i = 0; i < nparts; i++)
    {
      if(theta<=(i+1)*2*pi/nparts & theta > i*2*pi/nparts)
        nid[k] = i;
    }
  }
  //MF::mesh_node_partition(mesh, nparts, submeshes, nid, cid, "test_surface");
}

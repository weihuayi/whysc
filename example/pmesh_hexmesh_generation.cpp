#include <sstream> 
#include <cmath> 

#include <fstream>
#include <metis.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/HexJacobiQuality.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
typedef WHYSC::Mesh::ParallelMesh<GK, HexMesh> PMesh;
typedef WHYSC::Mesh::HexJacobiQuality<PMesh> MeshQuality;

typedef WHYSC::Mesh::VTKMeshReader<PMesh> Reader;
typedef WHYSC::Mesh::VTKMeshWriter<PMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;

const double pi = acos(-1.0);

int main(int argc, char * argv[])
{

  std::stringstream ss;
  ss << argv[1] << ".vtu";
  //auto mesh = std::make_shared<PQMesh>();
  auto mesh = std::make_shared<PMesh>();
  Reader reader(mesh);
  reader.read(ss.str());
  auto & gtag = mesh->get_node_int_data()["gtag"];
  auto & gdof = mesh->get_node_int_data()["gdof"];

  reader.get_node_data("gdof", gdof);
  reader.get_node_data("gtag", gtag);

  Writer writer(mesh);
  writer.set_points();
  writer.set_cells();
  writer.write("init_hex_mesh.vtu");

  auto & node = mesh->nodes();
  auto & cell = mesh->cells();
  for(auto & n : node)
    std::cout<< n <<std::endl;
  for(auto & c : cell)
    std::cout<< c[0] <<" " <<  c[1] << " " << c[2] << " " << c[3] << " " << c[4] <<" " <<  c[5] << " " << c[6] << " " << c[7]  <<std::endl;

  std::vector<PMesh> submeshes;

  Writer writer0(mesh);
  writer0.set_points();
  writer0.set_cells();
  writer0.write("init_mesh.vtu");

  int nparts = std::stoi(argv[2]);
  submeshes.resize(nparts);

  auto & nodes = mesh->nodes();
  auto NN = mesh->number_of_nodes();

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
  MF::mesh_node_partition(mesh, nparts, submeshes, nid, cid, "test_hex");
}

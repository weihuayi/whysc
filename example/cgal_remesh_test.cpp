#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <fstream>
#include <vector>
#include <sstream> 
#include <memory>
#include <map>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef Surface_mesh::Vertex_index vertex_descriptor;

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Toplogy Toplogy;
typedef WHYSC::Mesh::VTKMeshReader<TriMesh> Reader;
typedef WHYSC::Mesh::VTKMeshWriter Writer;
typedef WHYSC::Mesh::MeshFactory MF;

namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};

int read_file(const char* filename, Surface_mesh & mesh)
{
  std::ifstream input;
  input.open(filename);

	if (!input.is_open())
	{
    std::cout<< "读取文件失败" <<std::endl;
		return 1;
	}
  std::vector<vertex_descriptor> vdata;
  vdata.resize(1022);
  std::string data;

  int vsize = 1022;
  int csize = 2040;

  int i = 0;
  while(getline(input, data))
  {
    if(data[0] == 'f')
    {
      continue;
    }

    if(i < vsize+2) 
    {
      std::vector<double> v;
      std::string x="";
      for(int j = 0; j < data.length(); j++)
      {
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          v.push_back(stod(x));
          x="";
        }
      }
      vdata[i-1] = mesh.add_vertex(Point_3(v[0], v[1], v[2]));
      std::cout<< "herw " << mesh.point(vdata[i-1]) <<std::endl;
    } 
    else 
    {
      std::vector<int> c;
      std::string x="";
      for(int j = 0; j < data.length(); j++)
      {
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          c.push_back(stod(x));
          x="";
        }
      }
      std::cout<< "hahah" << c[0] << " " << c[1] << " " << c[2] <<std::endl;
      mesh.add_face(vdata[c[1]], vdata[c[2]], vdata[c[3]]);
    }
    i++;
  }
  return 0;
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  const char* filename1 = (argc > 1) ? argv[2] : "data/pig.off";
  Mesh mesh;

  /* 读文件并转换为 surface_mesh*/
  std::ifstream input;
  input.open(filename);

	if (!input.is_open())
	{
    std::cout<< "读取文件失败" <<std::endl;
		return 1;
	}
  std::vector<vertex_descriptor> vdata;
  vdata.resize(4532);
  std::string data;

  int vsize = 4532;

  int i = 0;
  while(getline(input, data))
  {
    if(data[0] == 'f')
    {
      i++;
      continue;
    }

    if(i <= vsize) 
    {
      std::vector<double> v;
      std::string x="";
      for(int j = 0; j < data.length(); j++)
      {
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          v.push_back(stod(x));
          x="";
        }
      }
      vdata[i-1] = mesh.add_vertex(Point_3(v[0], v[1], v[2]));
      std::cout<< i-1 << "herw " << mesh.point(vdata[i-1]) <<std::endl;
      i++;
    } 
    else 
    {
      std::vector<int> c;
      std::string x="";
      for(int j = 0; j < data.length(); j++)
      {
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          c.push_back(stoi(x)-1);
          x="";
        }
      }
      std::cout<< "hahah" << c[0] << " " << c[1] << " " << c[2] <<std::endl;
      mesh.add_face(vdata[c[0]], vdata[c[1]], vdata[c[2]]);
      i++;
    }
  }
  /*******************完成************************/
  TriMesh mesh0;
  std::cout<< "here" <<std::endl;
  MF::cgal_surface_mesh_to_triangle_mesh<Mesh, TriMesh>(mesh, mesh0);

  Writer writer;
  writer.set_points(mesh0);
  writer.set_cells(mesh0);
  writer.write("file.vtu");
 
  /*
  double target_edge_length = 0.04;
  unsigned int nb_iter = 3;
  std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh),
      mesh,
      boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
  std::cout << "done." << std::endl;
  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      .protect_constraints(true)//i.e. protect border, here
      );


  std::cout << "Remeshing done." << std::endl;

  std::ofstream out(filename1);
  out << "OFF"<< std::endl;

  std::map<int, int> indexmap;
  int kk = 1;
  for(auto & v : mesh.vertices())
  {
      auto i = v.idx();
      auto p = mesh.point(v);
      indexmap[i] = kk;
      std::cout<< " 编号 " << i <<std::endl;
      out << p.x() <<" " << p.y() << " "<< p.z() <<std::endl;
      kk++;
  }

  for( auto & f : mesh.faces())
  {
    auto i = f.idx();
    auto h0 = mesh.halfedge(f);
    auto h1 = mesh.next(h0);
    auto h2 = mesh.next(h1);

    out <<" " << indexmap[mesh.target(h0).idx()] <<" " << indexmap[mesh.target(h1).idx()] <<" " << indexmap[mesh.target(h2).idx()] <<std::endl;
    //out <<" " << mesh.target(h0).idx() <<" " << mesh.target(h1).idx() <<" " << mesh.target(h2).idx() <<std::endl;
  }
  out.close();
  */
} 

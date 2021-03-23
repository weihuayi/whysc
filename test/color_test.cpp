
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Mesh_triangulation_3.h>
//#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
//#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/Labeled_mesh_domain_3.h>
//#include <CGAL/make_mesh_3.h>
//#include <metis.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"
#include "mesh/MeshToplogy.h"
#include "mesh/color.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;
typedef WHYSC::Mesh::VTKMeshWriter<TriMesh> Writer;
typedef WHYSC::Mesh::MeshFactory MF;
typedef TriMesh::Toplogy Toplogy;


int main(int argc, char **argv)
{
   
    TriMesh mesh;
    MF::one_triangle_mesh(mesh);

    mesh.uniform_refine(7);
 
    int NN = mesh.number_of_nodes();
    int NC = mesh.number_of_cells();
    int NE = mesh.number_of_edges();

    int col[NN]={0};

    clock_t start, end;

    start = clock();
    mesh_coloring<TriMesh>(mesh, col);
    end = clock();

    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;

    //检验染色是否正确
    int a=0;
    Toplogy node2node;
    mesh.node_to_node(node2node);
    auto & nei = node2node.neighbors();
    auto & loc = node2node.locations();

    int cmin = 1000000000;
    int cmax = 0;
    for(int i = 0; i < NN; i++)
    {
      if(cmin>col[i])
      {
        cmin = col[i];
      }

      if(cmax<col[i])
      {
        cmax = col[i];
      }

      for(int j = loc[i]; j < loc[i+1]; j++){
        if(col[i]==col[nei[j]])
        {
          a++;
        }
      }
    }
    if(a==0)
    {
      std::cout<< "染色成功" << " " << "单元数" << " " << NC << " " << "节点数" << " " << NN << " " << "边数" << " " << NE <<endl; 
      std::cout<< "cmin " << cmin << " cmax " << cmax <<std::endl;
    }
    else
    {
      std::cout<< "染色失败" <<endl; 
    }

    //vtk网格顶点数据
    std::vector<double> nodedata;
    for(int i = 0; i < NN; i++)
    {
        nodedata.push_back((double)col[i]);
    }

    Writer writer(&mesh);
    writer.set_points();
    writer.set_point_data(nodedata, 1, "color");
    writer.set_cells();
    writer.write("color_test.vtu");

    return 0;
}

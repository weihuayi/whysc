#ifndef CGALMesher_2_h
#define CGALMesher_2_h

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

namespace WHYSC {
namespace Mesh {

class CGALMesher_2
{
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_2<K> Vb;
  typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
  typedef CDT::Vertex_handle Vertex_handle;
  typedef CDT::Point Point;
public:
  template<typename Box, typename TriMesh>
  static void boxmesh2d(Box & box, TriMesh & mesh)
  {
    CDT cdt;
  }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of CGALMesher_2_h

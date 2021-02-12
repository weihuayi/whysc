#ifndef MeshFactory_h
#define MeshFactory_h

#include <unordered_map>

namespace WHYSC {
namespace Mesh {

class MeshFactory
{

public:
    template<typename C3T3,  typename TetMesh>
    static void c3t3_to_tetmesh(C3T3 & c3t3, TetMesh & mesh)
    {
        typedef typename C3T3::Triangulation Tr;
        typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
        typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

        typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename Tr::Vertex_handle Vertex_handle;
        typedef typename Tr::Point Point; //can be weighted or not

        typedef typename TetMesh::Node Node;

        const Tr& tr = c3t3.triangulation();

        auto & nodes = mesh.nodes();
        nodes.resize(tr.number_of_vertices());

        std::unordered_map<Vertex_handle, int> V;
        int i = 0;
        for(auto vit = tr.finite_vertices_begin();
                vit != tr.finite_vertices_end();
                ++vit, ++i)
        {
            V[vit] = i;
            auto p = tr.point(vit);
            nodes[i][0] = p.x();
            nodes[i][1] = p.y();
            nodes[i][2] = p.z();
        }

        auto & cells = mesh.cells();
        cells.resize(c3t3.number_of_cells_in_complex());
        i = 0;
        for( auto cit = c3t3.cells_in_complex_begin();
           cit != c3t3.cells_in_complex_end();
           ++cit, ++i)
        {
            for(int j=0; j<4; j++)
              cells[i][j] = V[cit->vertex(j)];
        }
        mesh.init_top();
    }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshFactory_h

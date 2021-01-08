#ifndef Mesh_h
#define Mesh_h

#include <vector>
namespace WHYSC {

namespace Mesh {



template<typename GK, typename Vector, typename Node, typename Cell>
class Mesh
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;
    typedef typename Cell::Edge Edge;
    typedef typename Cell::Face Face;

    typedef std::vector<Node>::iterator NodeIterator;
    typedef std::vector<Edge>::iterator EdgeIterator;
    typedef std::vector<Face>::iterator FaceIterator;
    typedef std::vector<Cell>::iterator FaceIterator;

public:

    Mesh()
    {
        NN = 0;
        NE = 0;
        NF = 0;
        NC = 0;
    }

    I number_of_nodes()
    {
        return node.size();
    }

    I number_of_cells()
    {
        return cell.size();
    }

    I number_of_edges()
    {
        return edge.size();
    }

    I number_of_faces()
    {
        if(Cell::dim == 2)
            return edge.size();
        else if (Cell::dim == 3)
            return face.size();
    }

    int geo_dimension()
    {
        return Node::dim;
    }

    int top_dimension()
    {
        return Cell::dim;
    }

    void insert(Node & n)
    {
        node.push_back(n);
        NN += 1;
    }

    void insert(Cell & c)
    {
        cell.push_back(c);
        NC += 1;
    }

    void construct_top()
    {

    }



private:
    I NN;
    I NE;
    I NF;
    I NC;
    std::vector<Node> node;
    std::vector<Edge> edge;
    std::vector<Face> face;
    std::vector<Cell> cell;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Mesh_h

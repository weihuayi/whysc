#ifndef Mesh_h
#define Mesh_h

#include <vector>
#include <map>

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

    typedef typename std::vector<Node>::iterator NodeIterator;
    typedef typename std::vector<Edge>::iterator EdgeIterator;
    typedef typename std::vector<Face>::iterator FaceIterator;
    typedef typename std::vector<Cell>::iterator CellIterator;

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
        return Node::dimension();
    }

    int top_dimension()
    {
        return Cell::dimension();
    }

    void insert(Node n)
    {
        node.push_back(n);
        NN += 1;
    }

    void insert(Cell c)
    {
        cell.push_back(c);
        NC += 1;
    }

    void construct_top()
    {
        NF = 0;
        auto TD = Cell::dimension();
        std::map<I, I> idxmap;
        for(I i=0; i<NC; i++)
        {
            for(I j=0; j<Cell::ND[TD-1]; j++)
            {
               Face f;
               auto s = cell[i].get_local_face(j, f);
               auto it = idxmap.find(s);
               if(it == idxmap.end())
               {
                  cell[i].face[j] = NF;
                  idxmap.insert(std::pair<I, I>(s, NF));
                  f.face2cell[0] = i;
                  f.face2cell[1] = i;
                  f.face2cell[2] = j;
                  f.face2cell[3] = j;
                  face.push_back(f);
                  NF++;
               }
               else
               {
                  cell[i].face[j] = it->second;
                  face[it->second].face2cell[1] = i;
                  face[it->second].face2cell[3] = j;
               }
            }
        }

        if(TD == 3)
        {
            NE = 0;
            idxmap.clear();
            for(I i=0; i < NC; i++)
            {
                for(I j=0; j < Cell::ND[TD-2]; j++)
                {
                    Edge e;
                    auto s = cell[i].get_local_edge(j, e);
                    auto it = idxmap.find(s);
                    if(it == idxmap.end())
                    {
                        cell[i].edge[j] = NE;
                        idxmap.insert(std::pair<I, I>(s, NE));
                        edge.push_back(e);
                        NE++;
                    }
                    else
                    {
                        cell[i].edge[j] = it->second;
                    }
                }
            }
        }
    }

    void print()
    {
        for(I i = 0; i < NN; i++)
        {
            std::cout<< i << ":" <<  node[i] << std::endl;
        }

        for(I i = 0; i < NC; i++)
        {
            std::cout <<"cell"<<  i << ":" << cell[i][0] << " "<< cell[i][1] << " "
                << cell[i][2] << std::endl;
        }

        for(I i = 0; i < NF; i++)
        {
            std::cout <<"face"<<  i << ":" << face[i][0] << " "<< face[i][1] << " "
                << face[i][2] << std::endl;
        }
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

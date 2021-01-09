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

    typedef typename Cell::Face2Cell Face2Cell;
    typedef typename Cell::Cell2Face Cell2Face;

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
        return Node::dim;
    }

    int top_dimension()
    {
        return Cell::dim;
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
        face.resize(0);
        std::map<I, I> idxmap;
        I TD = Cell::dim;
        cell2face.resize(NC);
        I s = 0;
        //I localface[TD+1][TD] = Cell::face;

        for(I i=0; i<NC; i++)
        {
            Cell & c = cell[i];
            for(I j=0; j<TD+1; j++)
            {
               Face e;
               for(I k=0; k<TD; k++)
               {
                   e[k] = c[Cell::face[j][k]];
               }

               std::sort(e.begin(), e.end(), std::greater<I>());
               if(TD==2)
                   s = e[0] + e[1]*(e[1]+1)/2;
               if(TD==3)
                   s = e[0] + e[1]*(e[1]+1)/2 + e[2]*(e[2]+1)*(e[2]+2)/6;

               auto it = idxmap.find(s);
               if(it == idxmap.end())
               {
                  cell2face[i][j] = NF;
                  idxmap.insert(std::pair<I, I>(s, NF));
                  face.push_back(e);
                  Face2Cell f2c = {i, i, j, j};
                  face2cell.push_back(f2c);
                  NF++;
               }
               else
               {
                  cell2face[i][j] = it->second;
                  face2cell[it->second][1] = i;
                  face2cell[it->second][3] = j;
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
            std::cout << cell2face[i][0] << " " << cell2face[i][1] << " "
                << cell2face[i][2] << " " << std::endl;
        }

        for(I i = 0; i < NF; i++)
        {
            std::cout <<"face"<<  i << ":" << face[i][0] << " "<< face[i][1] << " "
                << face[i][2] << std::endl;
            std::cout << face2cell[i][0] << " " << face2cell[i][1] << " "
                << face2cell[i][2] << " " << face2cell[i][3] << std::endl;
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
    std::vector<Face2Cell> face2cell;
    std::vector<Cell2Face> cell2face;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Mesh_h

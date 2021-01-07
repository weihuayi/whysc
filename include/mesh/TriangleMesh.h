#ifndef TriangleMesh_2_h
#define TriangleMesh_2_h

#include <vector>
#include <map>

namespace WHYSC {
namespace Mesh {

template<typename GK>
class TriangleMesh 
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;

    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Vector_2 Vector_2;

    typedef typename std::array<I, 3> Cell;
    typedef typename std::array<I, 2> Edge;
    typedef typename std::array<I, 4> Edge2Cell;
    typedef typename std::array<I, 3> Cell2Edge;

    TriangleMesh()
    {
        NN = 0;
        NC = 0;
        NE = 0;
    }

    TriangleMesh(std::vector<F> nodes, std::vector<I> cells)
    {
       NN = nodes.size()/2;
       NC = cells.size()/3;
       node.resize(NN);
       cell.resize(NC);
       cell2edge.resize(NC);

       for(I i = 0; i < NN; i++)
       {
           node[i][0] = nodes[2*i];
           node[i][1] = nodes[2*i+1];
       }

       for(I i = 0; i < NC; i++)
       {
           cell[i][0] = cells[3*i];
           cell[i][1] = cells[3*i+1];
           cell[i][2] = cells[3*i+2];
       }
       construct_map();
    }

    void uniform_refine()
    {
       Cell oldCell[NC];
       for(I i=0; i<NC; i++)
       { 
           oldCell[i] = cell[i];
       }

       node.resize(NE+NN);
       cell.resize(NC*4);
       I n0, n1, n2, n3, n4, n5;
       for(I i=0; i<NE; i++)
       {
           n0 = edge[i][0];
           n1 = edge[i][1];

           node[NN+i] = 0.5*(node[n0]+node[n1]);
       }

       for(I i=0; i<NC; i++)
       {
           n0 = oldCell[i][0];
           n1 = oldCell[i][1];
           n2 = oldCell[i][2];
           n3 = cell2edge[i][0]+NN;
           n4 = cell2edge[i][1]+NN;
           n5 = cell2edge[i][2]+NN;

           cell[4*i] = {n3, n4, n5};
           cell[4*i+1] = {n0, n5, n4};
           cell[4*i+2] = {n1, n3, n5};
           cell[4*i+3] = {n2, n4, n3};
       }
       NN = NN+NE;
       NC = NC*4;
       construct_map();
    }

    void construct_map()
    {
        edge.resize(0);
        std::map<I, I> idxmap;
        I e0 = 0;
        I e1 = 0;
        I s = 0;
        NE = 0;
        cell2edge.resize(NC);
        I localEdge[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        for(I i=0; i < NC; i++)
        {
            Cell & c = cell[i];
            for(I j=0; j < 3; j++)
            {
                Edge e = {c[localEdge[j][0]], c[localEdge[j][1]]};
                if(e[0] < e[1])
                {
                    s = e[0] + e[1]*(e[1]+1)/2;
                }
                else
                {
                    s = e[1] + e[0]*(e[0]+1)/2;
                }
                auto it = idxmap.find(s);
                if(it == idxmap.end())
                {
                   cell2edge[i][j] = NE;
                   idxmap.insert(std::pair<I, I>(s, NE));
                   edge.push_back(e);
                   Edge2Cell e2c = {i, i, j, j};
                   edge2cell.push_back(e2c);
                   NE++;
                }
                else
                {
                   cell2edge[i][j] = it->second;
                   edge2cell[it->second][1] = i;
                   edge2cell[it->second][3] = j;
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

        for(I i = 0; i < NE; i++)
        {
            std::cout << "edge" << i << ":" << edge[i][0] << " " << edge[i][1] << " ";
            std::cout << "edge2cell"<< edge2cell[i][0] << " " << edge2cell[i][1] << " "
                << edge2cell[i][2] << " " << edge2cell[i][3] << std::endl;
        }

        for(I i = 0; i < NC; i++)
        {
            std::cout <<"cell"<<  i << ":" << cell[i][0] << " "<< cell[i][1] << " "
                << cell[i][2] << std::endl;
            std::cout << cell2edge[i][0] << " " << cell2edge[i][1] << " "
                << cell2edge[i][2] << " " << cell2edge[i][3] << std::endl;
        }
    }

private:
    I NN;
    I NE;
    I NC;
    std::vector<Point_2>  node;
    std::vector<Cell> cell;
    std::vector<Edge> edge;
    std::vector<Edge2Cell> edge2cell;
    std::vector<Cell2Edge> cell2edge;
};

} // end of namespace Mesh 

} // end of namespace WHYSC
#endif

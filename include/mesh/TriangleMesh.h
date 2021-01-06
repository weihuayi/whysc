#ifndef TriangleMesh_2_h
#define TriangleMesh_2_h

#include <vector>
#include <map>

namespace WHYSC {
namespace Mesh {

template<typename GK>
class TriangleMesh_2 
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

    Mesh()
    {
        NN = 0;
        NC = 0;
        NE = 0;
    }

    Mesh(std::vector<F> nodes, std::vector<I> cells)
    {
       NN = nodes.size()/2;
       NC = cells.size()/3;
       node.resize(NN);
       cell.resize(NC);

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
    }

    void construct_map()
    {
        std::map<I, I> idxmap;
        I e0 = 0;
        I e1 = 0;
        I s = 0;
        NE = 0;
        for(I i=0; i < NC; i++)
        {
            // 第 0 条边
            Cell & c = cell[i];
            e0 = c[1];
            e1 = c[2];
            if(e0 < e1)
            {
                 s = e0 + e1*(e1+1)/2;
            }
            else
            {
                s = e1 + e0*(e0+1)/2;
            }
            auto it = idxmap.find(s);
            if(it == idxmap.end())
            {
               cell2edge[i][0] = NE;
               idxmap.insert(std::pair<I, I>(s, NE));
               NE++;
            }
            else
            {
                cell2edge[i][0] = it->second;
            }


            // 第 1 条边
            e0 = c[2];
            e1 = c[0];
            if(e0 < e1)
            {
                s = e0 + e1*(e1+1)/2;
            }
            else
            {
                s = e1 + e0*(e0+1)/2;
            }
            auto it = idxmap.find(s);
            if(it == idxmap.end())
            {
               cell2edge[i][1] = NE;
               idxmap.insert(std::pair<I, I>(s, NE));
               NE++;
            }
            else
            {
                cell2edge[i][1] = it->second;
            }

            // 第 2 条边
            e0 = c[0];
            e1 = c[1];
            if(e0 < e1)
            {
                s = e0 + e1*(e1+1)/2;
            }
            else
            {
                s = e1 + e0*(e0+1)/2;
            }
            auto it = idxmap.find(s);
            if(it == idxmap.end())
            {
               cell2edge[i][2] = NE;
               idxmap.insert(std::pair<I, I>(s, NE));
               NE++;
            }
            else
            {
               cell2edge[i][2] = it->second;
            }
        }

        edge.resize(NE);
        edge2cell.resize(NE);
        for(I i=0; i < NC; i++)
        {
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

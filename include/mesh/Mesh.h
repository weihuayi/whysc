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
        m_NN = 0;
        m_NE = 0;
        m_NF = 0;
        m_NC = 0;
    }

    I number_of_nodes()
    {
        return m_nodes.size();
    }

    I number_of_cells()
    {
        return m_cells.size();
    }

    I number_of_edges()
    {
        return m_edges.size();
    }

    I number_of_faces()
    {
        if(Cell::dimension() == 2)
            return m_edges.size();
        else if (Cell::dimension() == 3)
            return m_faces.size();
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
        m_nodes.push_back(n);
        m_NN += 1;
    }

    void insert(Cell c)
    {
        m_cells.push_back(c);
        m_NC += 1;
    }

    void construct_top()
    {
        m_NF = 0;
        auto TD = Cell::dimension();
        std::map<I, I> idxmap;
        for(I i=0; i<m_NC; i++)
        {
            for(I j=0; j<Cell::ND[TD-1]; j++)
            {
               auto s = m_cells[i].local_face_index(j);
               auto it = idxmap.find(s);
               if(it == idxmap.end())
               {
                  m_cells[i].face(j) = m_NF;
                  idxmap.insert(std::pair<I, I>(s, m_NF));
                  m_faces.push_back(Face{i, i, j, j});
                  m_NF++;
               }
               else
               {
                  m_cells[i].face(j) = it->second;
                  m_faces[it->second].cell(1) = i;
                  m_faces[it->second].cell(3) = j;
               }
            }
        }

        if(TD == 3)
        {
            m_NE = 0;
            idxmap.clear();
            for(I i=0; i < m_NC; i++)
            {
                for(I j=0; j < Cell::ND[TD-2]; j++)
                {
                    auto s = m_cells[i].local_edge_index(j);
                    auto it = idxmap.find(s);
                    if(it == idxmap.end())
                    {
                        m_cells[i].edge(j) = m_NE;
                        idxmap.insert(std::pair<I, I>(s, m_NE));
                        m_edges.push_back(Edge{i, i, j, j});
                        m_NE++;
                    }
                    else
                    {
                        m_cells[i].edge(j) = it->second;
                    }
                }
            }
        }

        else if(TD == 2)
        {
            m_NE = m_NF;
        }
    }

    int * face(int & i)
    {
        int & N = Cell::NV[2];
        int *a = new int[N]; 

        Face & f = m_faces[i];
        Cell & c = m_cells[f.cell(0)];
        int & idx = f.index(0);

        for(int j=0; j<N; j++)
        {
            a[j] = c.node(Cell::localface[idx][j]);
        }
        return a;
    }

    int * edge(int & i)
    {
        int & N = Cell::NV[1];
        int *a = new int[N]; 

        Edge & e = m_edges[i];
        Cell & c = m_cells[e.cell(0)];
        int & idx = e.index(0);

        for(int j=0; j<N; j++)
        {
            a[j] = c.node(Cell::localedge[idx][j]);
        }
        return a;
    }


    void print()
    {
        for(I i = 0; i < m_NN; i++)
        {
            std::cout<< "point " << i << ":" <<  m_nodes[i] << std::endl;
        }

        for(I i = 0; i < m_NF; i++)
        {
            std::cout<< "face2node " << i << ":";
            for(int j = 0; j < Cell::NV[2]; j++)
            {
                std::cout<<  " " << face(i)[j];
            }
            std::cout << std::endl;
        }

        for(I i = 0; i < m_NE; i++)
        {
            std::cout<< "edge2node " << i << ":";
            for(int j = 0; j < Cell::NV[1]; j++)
            {
                std::cout<<  " " << edge(i)[j];
            }
            std::cout << std::endl;
        }

        auto TD = Cell::dimension();
        for(I i = 0; i < m_NC; i++)
        {
            std::cout <<"cell2node "<<  i << ":";
            for(I j = 0; j < Cell::ND[0]; j++)
                std::cout <<  " " << m_cells[i].node(j);
            std::cout << std::endl;
        }

        for(I i = 0; i < m_NC; i++)
        {
            std::cout <<"cell2edge "<<  i << ":";
            for(I j = 0; j < Cell::ND[1]; j++)
                std::cout << " " << m_cells[i].edge(j);
            std::cout << std::endl;
        }

        for(I i = 0; i < m_NC; i++)
        {
            std::cout <<"cell2face "<<  i << ":";
            for(I j = 0; j < Cell::ND[2]; j++)
                std::cout << " " << m_cells[i].face(j);
            std::cout << std::endl;
        }
    }


private:
    I m_NN;
    I m_NE;
    I m_NF;
    I m_NC;
    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    std::vector<Face> m_faces;
    std::vector<Cell> m_cells;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Mesh_h

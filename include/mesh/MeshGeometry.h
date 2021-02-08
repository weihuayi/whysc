#ifndef MeshGeometry_h
#define MeshGeometry_h


#include <vector>
#include <cassert>

namespace WHYSC {
namespace Mesh {

/*
 * MeshGeometry:  网格几何对象, 负责几何解释浮点数组中存储的节点坐标
 *
 */
template<typename GK>
class MeshGeometry
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;

    typedef typename GK::PPoint_3 Node;
    typedef typename GK::Point_3 Point;
    typedef typename GK::Vector_3 Vector;

    class NodeIterator
    {
    public:
        NodeIterator(F * data, const I NN, const I id=0)
        {
            m_data = data;
            m_id = id;
            m_NN = NN;
        }

        ~NodeIterator(){};

        bool end() {return m_id == m_NN;}

        void operator++() // prefix
        {
            assert(m_id < m_NN);
            m_data += 3;
            ++m_id;
        }
        
        void operator++(int) // postfix
        {
            assert(m_id < m_NN);
            m_data += 3;
            ++m_id;
        }

        void operator--() // prefix
        {
            assert(m_id > 0);
            m_data -= 3;
            --m_id;
        }
        
        void operator--(int) // postfix
        {
            assert(m_id > 0);
            m_data -= 3;
            m_id--;
        }

        Node value() { return Node(m_data);}

        Node operator*()
        {
            return Node(m_data);
        }

        bool operator!=(const NodeIterator & it)
        {
            return m_id != it.m_id;
        }

        bool operator==(const NodeIterator & it)
        {
            return m_id == it.m_id;
        }

    private:
        F * m_data;
        I m_id;
        I m_NN;
    };

public:
    MeshGeometry()
    {
        m_node = NULL;
        m_NN = 0;
        m_NR = 0;
        m_owner = true;
    }

    MeshGeometry(I NR)
    {
        m_node = new F[3*NR];
        m_owner = true;
        m_NN = 0;
        m_NR = NR;
    }

    MeshGeometry(const F * node, const I NN, bool owner=false)
    {
        m_owner = owner;
        if(m_owner)
        {
            m_node = new F[3*NN]; 
            m_NN = NN;
            m_NR = NN;
        }
        else
        {
            m_node = node;
            m_NN = NN;
            m_NR = NN;
        }
    }


    void clear()
    {
        if(m_owner && m_node != NULL)
        {
            delete [] m_node;
        }
        m_node = NULL;
        m_NN = 0;
        m_NR = 0;
        m_owner = true;
    }

    ~MeshGeometry()
    {
        if(m_owner && m_node != NULL)
            delete [] m_node;
    }

    F * data()
    {
        return m_node;
    }

    void reseize(I NN)
    {
    }

    void reserve(I NR)
    {
    }

    /*
     *
     * Notes
     * -----
     *  在第 i 个位置插入一个节点的坐标 
     *
     */
    void insert(const I i, const F x, const F y, const F z=0.0)
    {
        if(i < m_NN)
        {
            m_node[3*i] = x;
            m_node[3*i + 1] = y;
            m_node[3*i + 2] = z;
        }
        else
        {
            //TODO: resize memory 
        }
    }

    /*
     *
     * Notes
     * -----
     *  在尾部增加一个节点的坐标
     *
     */
    void push_back(const F x, const F y, const F z=0.0)
    {
        if(m_NN < m_NR)
        {
            m_node[3*m_NN] = x;
            m_node[3*m_NN + 1] = y;
            m_node[3*m_NN + 2] = z;
            m_NN += 1;
        }
        else
        {
            //TODO: resize memory 
        }
    }

    NodeIterator node_begin()
    {
        return NodeIterator(m_node, m_NN);
    }

    NodeIterator node_end()
    {
        return NodeIterator(m_node, m_NN, m_NN);
    }

    I number_of_nodes()
    {
        return m_NN;
    }

    static int geo_dimension()
    {
        return 3;
    }

    Node & operator [] (const int i)
    {
        return Node(m_node + 3*i);
    }

    Node node(I i)
    {
        return Node(m_node + 3*i);
    }

    void node(std::vector<Node> & nodes)
    {
        nodes.resize(m_NN);
        for(I i = 0; i < m_NN; i++)
            nodes[i] = m_node + 3*i;
    }

    void print()
    {
        std::cout << "The total number of nodes is " << m_NN << std::endl;
        for(I i = 0; i < m_NN; i++)
        {
            auto n = node(i);
            std::cout << i << ": " << n << std::endl;
        }
    }

private:
    bool m_owner;
    F * m_node;
    I m_NN; // 实际存储的节点个数 
    I m_NR; // 实际分配内存的个数, 在任何情况下, m_NN <= m_NR
};
} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of MeshGeometry_h

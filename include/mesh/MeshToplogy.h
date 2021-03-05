#ifndef MeshToplogy_h
#define MeshToplogy_h

namespace WHYSC {
namespace Mesh {
/*
 * MeshToplogy: 定义一个网格中实体 A 和 实体 B 之间的拓扑关系
 *
 */
template<typename I, typename Container=std::vector<I> >
class MeshToplogy
{
public:
  typedef typename Container::reference         Reference;
  typedef typename Container::const_reference   ConstReference;
  typedef typename Container::iterator          Iterator;
  typedef typename Container::const_iterator    ConstIterator;

  class AdjEntitiesIterator
  {
  public:

    AdjEntitiesIterator(
        Iterator offset_it, 
        Iterator adj_begin,
        Iterator loc_begin):
      m_offset_it(offset_it), 
      m_adj_begin(adj_begin),
      m_loc_begin(loc_begin) {}

    ~AdjEntitiesIterator(){}

    void operator++() // prefix
    {
      ++m_offset_it;
    }

    void operator++(int) //postfix
    {
      m_offset_it++;
    }

    size_t size()
    {
      return *(m_offset_it+1) - *(m_offset_it);
    }

    Reference operator[](std::size_t i)
    {
      return m_adj_begin[*m_offset_it + i];
    }

    Iterator adj_begin()
    {
      return m_adj_begin[*m_offset_it];
    }

    Iterator adj_end()
    {
      return m_adj_begin[*(m_offset_it + 1)];
    }

    Iterator loc_begin()
    {
      return m_loc_begin[*m_offset_it];
    }

    Iterator loc_end()
    {
      return m_loc_begin[*(m_offset_it + 1)];
    }

    bool operator==(const AdjEntitiesIterator& it) const
    {
      return (m_offset_it == it.m_offset_it) && (m_adj_begin == m_adj_begin); 
    }

    bool operator!=(const AdjEntitiesIterator& it) const
    {
      return (m_offset_it != it.m_offset_it) || (m_adj_begin != m_adj_begin); 
    }

  private:
    Iterator m_offset_it;
    Iterator m_adj_begin;
    Iterator m_loc_begin;
  };

public:
  MeshToplogy()
  {
    m_TA = -1;
    m_TB = -1;
    m_NA = 0;
    m_NB = 0;
  }

  MeshToplogy(I TA, I TB, I NA, I NB)
  {
    init(TA, TB, NA, NB);
  }

  void init(I TA, I TB, I NA, I NB)
  {
    m_TA = TA;
    m_TB = TB;
    m_NA = NA;
    m_NB = NB;
    m_offset.resize(m_NA+1);
  }

  I & top_dimension_of_entity_A()
  {
    return m_TA;
  }

  I & top_dimension_of_entity_B()
  {
    return m_TB;
  }

  bool empty()
  {
    return (m_TA == -1) || (m_TB == -1) || (m_NA == 0) || (m_NB == 0);
  }

  void clear()
  {
    m_adj.clear();
    m_loc.clear();
    m_offset.clear();
    m_TA = -1;
    m_TB = -1;
    m_NA = 0;
    m_NB = 0;
  }

public: // new interface

    I number_of_adj_entities(const I i)
    {
      return m_offset[i+1] - m_offset[i];
    }

    I * adj_entities(const I i)
    {
      return &m_adj[m_offset[i]];
    }

    AdjEntitiesIterator operator[](std::size_t i)
    {
      return AdjEntitiesIterator(m_offset.begin()+i, m_adj.begin(), m_loc.begin());
    }

    AdjEntitiesIterator adj_begin()
    {
      return AdjEntitiesIterator(m_offset.begin(), m_adj.begin(), m_loc.begin());
    }

    AdjEntitiesIterator adj_end()
    {
      return AdjEntitiesIterator(--m_offset.end(), m_adj.begin(), m_loc.begin());
    }

public: // deprecated interface

    Container & neighbors()
    {
        return m_adj;
    }

    Container & local_indices()
    {
        return m_loc;
    }

    Container & locations()
    {
        return m_offset;
    }

    I * neighbors(const I i)
    {
      return &m_adj[m_offset[i]];
    }

    I number_of_neighbors(const I i)
    {
      return m_offset[i+1] - m_offset[i];
    }


private:
    Container m_adj; // 存储每个 A 实体相邻的 B 实体编号
    Container m_loc; // 存储每个 A 实体在相邻的 B 实体中的局部编号
    Container m_offset; // 每个 A 实体邻接实体的偏移量 
    I m_TA; // A 实体的拓扑维数
    I m_TB; // B 实体的拓扑维数
    I m_NA; // the number of A entity
    I m_NB; // the number of B entity
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshToplogy_h

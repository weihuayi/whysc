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

  class AdjEntitySet 
  {
    public:
      AdjEntitySet(I num, I * adj)
      {
        m_num = num;
        m_adj = adj;
      }

      I & operator [] (const int i)
      {
          return m_adj[i];
      }

      const I & operator [] (const int i) const
      {
          return m_adj[i];
      }

      I adj_entity(const int i)
      {
        return m_adj[i];
      }
      
    private:
      I m_num;
      I * m_adj;
  };

  class AdjEntitySetWithLoc
  {
    public:
      AdjEntitySetWithLoc(I num, I * adj, I * loc)
      {
        m_num = num;
        m_adj = adj;
        m_loc = loc;
      }

      I adj_entity(const int i)
      {
        return m_adj[i];
      }

      I adj_local_index(const int i)
      {
        return m_loc[i];
      }

    private:
      I m_num;
      I * m_adj;
      I * m_loc;
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

    /*  第 i 个 A 实体相邻的 B 实体的个数 */
    I number_of_adj_entities(const I i)
    {
      return m_offset[i+1] - m_offset[i];
    }

    I * adj_entities_begin(const I i)
    {
      return &m_adj[m_offset[i]];
    }

    I * adj_local_index(const I i)
    {
      return &m_loc[m_offset[i]];
    }

    AdjEntitySet adj_entities(const I i)
    {
      return AdjEntitySet(i, &m_adj[m_offset[i]]);
    }

    AdjEntitySetWithLoc adj_entities_with_local(const I i)
    {
      return AdjEntitySetWithLoc(i, &m_adj[m_offset[i]], &m_loc[m_offset[i]]);
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
    Container m_loc; // 存储每个 A 实体在相邻的 B 实体中的局部编号, 如 A 是拓扑维数高于 B 的实体， m_loc 是空
    Container m_offset; // 每个 A 实体邻接实体的偏移量 
    I m_TA; // A 实体的拓扑维数
    I m_TB; // B 实体的拓扑维数
    I m_NA; // the number of A entity
    I m_NB; // the number of B entity
};



} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshToplogy_h

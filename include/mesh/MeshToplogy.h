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
    MeshToplogy()
    {
        m_NA = 0;
        m_NB = 0;
    }

    MeshToplogy(I na, I nb)
    {
        m_NA = na;
        m_NB = nb;
        m_location.resize(m_NA+1);
    }
private:
    Container m_neighbor;
    Container m_localidx;
    Container m_location;
    I m_NA; // the number of A entity
    I m_NB; // the number of B entity
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshToplogy_h

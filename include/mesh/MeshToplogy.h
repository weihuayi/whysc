#ifndef Toplogy_h
#define Toplogy_h

namespace WHYSC {
namespace Mesh {
/*
 * MeshToplogy: 定义一个网格中实体 A 和 实体 B 之间的拓扑关系
 *
 */
template<typename I, typename Container=std::vector<I> >
class Toplogy
{
public:
    MeshToplogy()
    {
        m_NA = 0;
        m_NB = 0;
    }

    MeshToplogy(I NA, I NB)
    {
        m_NA = NA;
        m_NB = NB;
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
#endif // end of Toplogy_h

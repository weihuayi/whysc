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
        m_TA = -1;
        m_TB = -1;
        m_NA = 0;
        m_NB = 0;
    }

    MeshToplogy(I TA, I TB, I NA, I NB)
    {
        m_TA = TA;
        m_TB = TB;
        m_NA = NA;
        m_NB = NB;
        m_location.resize(m_NA+1);
    }

private:
    Container m_neighbor; // 存储每个 A 实体相邻的 B 实体编号
    Container m_localidx; // 存储每个 A 实体在相邻的 B 实体中的局部编号 
    Container m_location; //
    I m_TA; // A 实体的拓扑维数
    I m_TB; // B 实体的拓扑维数
    I m_NA; // the number of A entity
    I m_NB; // the number of B entity
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of MeshToplogy_h

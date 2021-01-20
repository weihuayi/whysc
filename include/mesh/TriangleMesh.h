#ifndef TriangleMesh_h
#define TriangleMesh_h

#include <vector>
#include <map>

namespace WHYSC {
namespace Mesh {

template<typename GK, typename Node, typename Vector, typename Container=std::vector>
class TriangleMesh 
{
public:
    typedef typename GK::Int I;
    typedef typename GK::Float F;

    TriangleMesh()
    {
    }

private:
    Container<Node>  m_nodes;
    Container<Cell> cell;
};

} // end of namespace Mesh 

} // end of namespace WHYSC

#endif // end of TriangleMesh_h

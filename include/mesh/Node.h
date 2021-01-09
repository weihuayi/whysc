#ifndef Node_h
#define Node_h

namespace WHYSC {

namespace Mesh {

template<typename Point>
class Node: public Point
{
private:
    int m_id;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Node_h

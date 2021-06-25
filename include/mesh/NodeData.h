#ifndef NodeData_h
#define NodeData_h

#include <string>
#include <map>
#include <vector>

namespace WHYSC {
namespace Mesh {


template<typename T>
class NodeData: public std::map<std::string, std::vector<T> >
{
};

typedef NodeData<int> NodeIntData;
typedef NodeData<double> NodeDoubleData;

// NodeIntData["gdof"].resize()

}
}

#endif // end of NodeData_h

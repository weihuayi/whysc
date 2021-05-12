#ifndef TriangleMeshQuality_h
#define TriangleMeshQuality_h

#include <vector>
#include <array>
#include <map>

#include "TriangleMesh.h"

using json = nlohmann::json;

namespace WHYSC {
namespace Mesh {

class TriangleMeshQuality
{
public:
  TriangleMeshQuality(std::shared_ptr<TriangleMesh> mesh)
  {
    m_mesh = mesh;
  }

  void norm_quality(std::vector<double> & quality)
  {
    auto NC = m_mesh->number_of_cells();
    auto & cell = m_mesh->cells();
    auto & node = m_mesh->nodes();

    quality.resize(NC);
    for(int i = 0; i < NC; i++)
    {
      auto v0 = node[cell[i][1]] - node[cell[i][0]];
      auto v1 = node[cell[i][2]] - node[cell[i][1]];
      auto v2 = node[cell[i][0]] - node[cell[i][2]];

      auto L0 = v0.squared_length(); 
      auto L1 = v1.squared_length();
      auto L2 = v2.squared_length();

      auto S = std::sqrt(cross(v0, v1).squared_length());
      auto r = S/(std::sqrt(L0)+std::sqrt(L1)+std::sqrt(L2));

      auto C = std::pow(dot(v0, v2), 2)/L0*L2;
      auto R = std::sqrt(L1/(4*(1-C)));
      quality[i] = r/R;
    }
  }

private:
  std::shared_ptr<TriangleMesh> m_mesh;
};
} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of TriangleMeshQuality_h

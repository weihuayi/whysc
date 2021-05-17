#ifndef GObject_h
#define GObject_h
/*
 * 我们想设计一个描述几何区域的类， 一个几何区域由三种基本几何实体组成
 * - 顶点: 0 维几何实体
 * - 曲线: 1 维几何实体
 * - 曲面: 2 维几何实体
 * - 部件: 3 维几何实体
 */
namespace WHYSC {
namespace GeometryObject {
    
enum GObjectType 
{
    VERTEX = 0,
    CURVE = 1,
    SURFACE = 2,
    PART = 3
};

class GObject 
{
public:
  GObject(const std::string & name, const int id, ): m_name(name), m_id(id)
  {
  }
  const std::string& name() const { return m_name; }
  const int id() const {return m_id;}
  const GObjectType type() const {return m_type;}
protected:
    std::string m_name;
    int m_id;
    GObjectType m_type;
};
} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of GObject_h
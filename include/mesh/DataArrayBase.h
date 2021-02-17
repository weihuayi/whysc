#ifndef DataArrayBase_h
#define DataArrayBase_h

#include <string>
#include <typeinfo>

namespace WHYSC {
namespace Mesh {

class DataArrayBase
{
public:
  DataArrayBase(const std::string& name): m_name(name) {}
public:
  virtual ~DataArrayBase() {};
  virtual size_t size() = 0;
  virtual void reserve(size_t n) = 0;
  virtual void resize(size_t n) = 0;
  virtual void shrink_to_fit() = 0;
  virtual void push_back() = 0;
  virtual void reset(size_t idx) = 0;
  virtual bool transfer(const DataArrayBase& other) = 0;
  virtual bool transfer(const DataArrayBase& other, std::size_t from, std::size_t to) = 0;
  virtual void swap(size_t i0, size_t i1) = 0;
  virtual DataArrayBase* clone () const = 0;
  virtual DataArrayBase* empty_clone () const = 0;
  virtual const std::type_info& type() const = 0;
public:
  const std::string& name() const { return m_name; }
  bool is_same (const DataArrayBase& other)
  {
    return (name() == other.name() && type() == other.type());
  }
protected:
  std::string m_name;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of DataArrayBase_h

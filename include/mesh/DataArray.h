#ifndef DataArray_h
#define DataArray_h

#include <vector>
#include "DataArrayBase.h"

namespace WHYSC {
namespace Mesh {

template<typename T>
class DataArray: DataArrayBase
{
public:
  typedef T ValueType;
  typedef std::vector<ValueType> VectorType;
  typedef typename VectorType::reference         Reference;
  typedef typename VectorType::const_reference   ConstReference;
  typedef typename VectorType::iterator          Iterator;
  typedef typename VectorType::const_iterator    ConstIterator;

  DataArray(const std::string& name, T t=T()): DataArrayBase(name), m_value(t) {}
public:

  virtual size_t size()
  {
    return m_data.size();
  }

  virtual void reserve(size_t n)
  {
    m_data.reserve(n);
  }

  virtual void resize(size_t n)
  {
    m_data.resize(n, m_value);
  }

  virtual void push_back()
  {
    m_data.push_back(m_value);
  }

  virtual void reset(size_t idx)
  {
    m_data[idx] = m_value;
  }

  bool transfer(const DataArrayBase& other)
  {
    const DataArray<T>* p = dynamic_cast<const DataArray<T>*>(&other);
    if(p != nullptr){
      std::copy((*p).m_data.begin(), (*p).m_data.end(), m_data.end()-(*p).m_data.size());
      return true;
    }
    return false;
  }

  bool transfer(const DataArrayBase& other, std::size_t from, std::size_t to)
  {
    const DataArray<T>* p = dynamic_cast<const DataArray<T>*>(&other);
    if (p != nullptr)
    {
      m_data[to] = (*p)[from];
      return true;
    }
    return false;
  }

  virtual void shrink_to_fit()
  {
    VectorType(m_data).swap(m_data);
  }

  virtual void swap(size_t i0, size_t i1)
  {
    T d(m_data[i0]);
    m_data[i0] = m_data[i1];
    m_data[i1] = d;
  }

  virtual DataArrayBase * clone() const
  {
    DataArray<T>* p = new DataArray<T>(this->m_name, this->m_value); 
    p->m_data = m_data;
    return p;
  }

  virtual DataArrayBase * empty_clone() const
  {
    DataArray<T>* p = new DataArray<T>(this->m_name, this->m_value); 
    return p;
  }

  virtual const std::type_info& type() const { return typeid(T);}

public:

  const T* data() const
  { // does not work for T==bool
    return &m_data[0];
  }

  Reference operator[](std::size_t i)
  {
    return m_data[i];
  }

  ConstReference operator[](std::size_t i) const
  {
    return m_data[i];
  }

  Iterator begin() { return m_data.begin();}
  Iterator end() {return m_data.end();}
  ConstIterator begin() const { return m_data.begin();}
  ConstIterator end() const {return m_data.end();}

private:
  VectorType m_data;
  ValueType m_value;
};

typedef DataArray<double> DoubleDataArray;
typedef DataArray<int> IntDataArray;

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of DataArray_h

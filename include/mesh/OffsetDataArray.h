#ifndef OffsetDataArray_h
#define OffsetDataArray_h

#include <vector>
#include "DataArrayBase.h"

namespace WHYSC {
namespace Mesh {

template<typename T, typename I=int>
class OffsetDataArray: DataArrayBase
{
public:

  typedef I Index;
  typedef T ValueType;
  typedef std::vector<ValueType> VectorType;
  typedef typename VectorType::reference         Reference;
  typedef typename VectorType::const_reference   ConstReference;
  typedef typename VectorType::iterator          Iterator;
  typedef typename VectorType::const_iterator    ConstIterator;

  typedef std::vector<Index> OffsetType;

  class OffsetIterator
  {
  public:
    OffsetIterator(VectorType * data, OffsetType * offset):m_begin(begin), m_end(end), m_off(dataend) {}
    ~OffsetIterator(){}

    void operator++() // prefix
    {
      m_begin += m_csize;
    }

    void operator++(int) //postfix
    {
      m_begin += m_csize;
    }

    Reference operator[](std::size_t i)
    {// 0 <= i < m_csize
      return *(m_begin + i);
    }

    bool end() {return m_begin == m_end;}

    bool operator==(const ComponentIterator& it)
    {
      return (m_begin == it.m_begin) && (m_end == it.m_end) && (m_csize == m_csize); 
    }

    bool operator!=(const OffsetIterator& it)
    {
      return (m_begin != it.m_begin) || (m_end != it.m_end) || (m_csize != m_csize);
    }
  private:
    Iterator m_begin;
    OffsetType* m_offset;
  };

  ComponentDataArray(const std::string& name, int csize=1, T t=T()): DataArrayBase(name), m_value(t), m_csize(csize) {}

public:

  virtual size_t size()
  {
    return m_data.size();
  }

  virtual void reserve(size_t n)
  {
    m_data.reserve(n*m_csize);
  }

  virtual void resize(size_t n)
  {
    m_data.resize(n*m_csize, m_value);
  }

  virtual void push_back()
  {
    for(int i = 0; i < m_csize; i++)
      m_data.push_back(m_value);
  }

  virtual void reset(size_t idx)
  {
    for(int i = 0; i < m_csize; i++)
      m_data[idx*m_csize + i] = m_value;
  }

  bool transfer(const DataArrayBase& other)
  {
    const ComponentDataArray<T>* p = dynamic_cast<const ComponentDataArray<T>*>(&other);
    if(p != nullptr){
      std::copy((*p).m_data.begin(), (*p).m_data.end(), m_data.end()-(*p).m_data.size());
      return true;
    }
    return false;
  }

  bool transfer(const DataArrayBase& other, std::size_t from, std::size_t to)
  {
    const ComponentDataArray<T>* p = dynamic_cast<const ComponentDataArray<T>*>(&other);
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
    for(int i=0; i < m_csize; i++)
    {
      T d(m_data[i0*m_csize + i]);
      m_data[i0*m_csize + i] = m_data[i1*m_csize + i];
      m_data[i1*m_csize + i] = d;
    }
  }

  virtual DataArrayBase * clone() const
  {
    ComponentDataArray<T>* p = new ComponentDataArray<T>(this->m_name, this->m_size, this->m_value); 
    p->m_data = m_data;
    return p;
  }

  virtual DataArrayBase * empty_clone() const
  {
    DataArray<T>* p = new DataArray<T>(this->m_name, this->m_size, this->m_value); 
    return p;
  }

  virtual const std::type_info& type() const { return typeid(T);}

public:

  size_t number_of_components()
  {
    return m_data.size()/m_csize;
  }

  int component_size() { return m_csize;}

  ComponentIterator component_begin()
  {
    return ComponentIterator(m_data.begin(), m_data.end(), m_csize);
  }

  ComponentIterator component_end()
  {
    return ComponentIterator(m_data.end(), m_data.end(), m_csize);
  }

  Iterator component_begin(std::size_t i) { return m_data.begin()+i*m_csize;}
  Iterator component_end(std::size_t i) {return m_data.begin()+(i+1)*m_csize;}
  ConstIterator component_begin(std::size_t i) const { return m_data.begin()+i*m_csize;}
  ConstIterator component_end(std::size_t i) const {return m_data.begin()+(i+1)*m_csize;}

private:
  VectorType m_data;
  ValueType m_value;
  int m_csize; // the size of each component
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of OffsetDataArray_h

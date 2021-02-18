#ifndef ComponentDataArray_h
#define ComponentDataArray_h

#include <vector>
#include "DataArrayBase.h"

namespace WHYSC {
namespace Mesh {

template<typename T>
class ComponentDataArray: DataArrayBase
{
public:

  typedef T ValueType;
  typedef std::vector<ValueType> VectorType;
  typedef typename VectorType::reference         Reference;
  typedef typename VectorType::const_reference   ConstReference;
  typedef typename VectorType::iterator          Iterator;
  typedef typename VectorType::const_iterator    ConstIterator;

  class ComponentIterator
  {
  public:
    ComponentIterator(Iterator begin, Iterator end, int csize):m_begin(begin), m_end(end), m_csize(end) {}
    ~ComponentIterator(){}

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

    bool operator==(const ComponentIterator& it) const
    {
      return (m_begin == it.m_begin) && (m_end == it.m_end) && (m_csize == m_csize); 
    }

    bool operator!=(const ComponentIterator& it) const
    {
      return (m_begin != it.m_begin) || (m_end != it.m_end) || (m_csize != m_csize);
    }
  private:
    Iterator m_begin;
    Iterator m_end;
    int m_csize;
  };

  ComponentDataArray(const std::string& name, int csize=1, T t=T()): DataArrayBase(name), m_value(t), m_csize(csize) {}

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
    T d(m_data[i0]);
    m_data[i0] = m_data[i1];
    m_data[i1] = d;
  }

  virtual DataArrayBase * clone() const
  {
    ComponentDataArray<T>* p = new ComponentDataArray<T>(this->m_name, this->m_csize, this->m_value); 
    p->m_data = m_data;
    return p;
  }

  virtual DataArrayBase * empty_clone() const
  {
    ComponentDataArray<T>* p = new ComponentDataArray<T>(this->m_name, this->m_csize, this->m_value); 
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

public:

  size_t number_of_components()
  {
    return m_data.size()/m_csize;
  }

  int component_size(size_t i) { return m_csize;}

  ComponentIterator component_begin()
  {
    return ComponentIterator(m_data.begin(), m_data.end(), m_csize);
  }

  ComponentIterator component_end()
  {
    return ComponentIterator(m_data.end(), m_data.end(), m_csize);
  }

  Iterator component_begin(size_t i) { return m_data.begin()+i*m_csize;}
  Iterator component_end(size_t i) {return m_data.begin()+(i+1)*m_csize;}
  ConstIterator component_begin(size_t i) const { return m_data.begin()+i*m_csize;}
  ConstIterator component_end(size_t i) const {return m_data.begin()+(i+1)*m_csize;}

private:
  VectorType m_data;
  ValueType m_value;
  int m_csize; // the size of each component
};


} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of ComponentDataArray_h

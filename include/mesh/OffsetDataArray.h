#ifndef OffsetDataArray_h
#define OffsetDataArray_h

#include <vector>
#include "DataArray.h"

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
  typedef typename OffsetType::iterator OffsetIterator;
  typedef typename OffsetType::const_iterator ConstOffsetIterator;

  class ComponentIterator
  {
  public://注意这里的 offset_end 要指向 offset 数组的最后一个元素
    ComponentIterator(OffsetIterator offset_it, OffsetIterator offset_last, Iterator data_begin):
      m_offset_it(offset_it), m_offset_last(offset_last), m_data_begin(data_begin) {}
    ~ComponentIterator(){}

    void operator++() // prefix
    {
      ++m_offset_it;
    }

    void operator++(int) //postfix
    {
      m_offset_it++;
    }

    bool end() {return m_offset_it == m_offset_last;}

    size_t size()
    {
      return *(m_offset_it+1) - *(m_offset_it);
    }

    Reference operator[](std::size_t i)
    {
      return m_data_begin[*m_offset_it + i];
    }

    Iterator component_begin()
    {
      return m_data_begin[*m_offset_it];
    }

    Iterator component_end()
    {
      return m_data_begin[*(m_offset_it + 1)];
    }

    bool operator==(const ComponentIterator& it) const
    {
      return (m_offset_it == it.m_offset_it) && (m_offset_last == it.m_offset_last) && (m_data_begin == m_data_begin); 
    }

    bool operator!=(const ComponentIterator& it) const
    {
      return (m_offset_it != it.m_offset_it) || (m_offset_last != it.m_offset_last) || (m_data_begin != m_data_begin); 
    }

  private:
    OffsetIterator m_offset_it;
    OffsetIterator m_offset_last;
    Iterator m_data_begin;
  };

  OffsetDataArray(const std::string& name, T t=T()): DataArrayBase(name), m_value(t) {}

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
    const OffsetDataArray<T>* p = dynamic_cast<const OffsetDataArray<T>*>(&other);
    if(p != nullptr)
    {
      std::copy((*p).m_data.begin(), (*p).m_data.end(), m_data.end()-(*p).m_data.size());
      return true;
    }
    return false;
  }

  bool transfer(const DataArrayBase& other, std::size_t from, std::size_t to)
  {
    const OffsetDataArray<T>* p = dynamic_cast<const OffsetDataArray<T>*>(&other);
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
    OffsetDataArray<T>* p = new OffsetDataArray<T>(this->m_name, this->m_value); 
    p->m_data = m_data;
    p->m_offset = m_offset;
    return p;
  }

  virtual DataArrayBase * empty_clone() const
  {
    OffsetDataArray<T>* p = new OffsetDataArray<T>(this->m_name, this->m_value); 
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

  OffsetType& offset() 
  {
    return m_offset;
  }

  size_t number_of_components()
  {
    return m_offset.size()-1;
  }

  ComponentIterator component_begin()
  {
    return ComponentIterator(m_offset.begin(), m_offset.end()--, m_data.begin());
  }

  ComponentIterator component_end()
  {
    return ComponentIterator(m_offset.begin(), m_offset.end()--, m_data.begin());
  }

private:
  VectorType m_data;
  OffsetType m_offset;
  ValueType m_value;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of OffsetDataArray_h

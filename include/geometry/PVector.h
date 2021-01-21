#ifndef PVector_h
#define PVector_h

namespace WHYSC {
namespace GeometryObject {

template<typename F, int DIM>
class PVector
{
public:
    typedef F Float;
public:

    PVector(F * p)
    {
        m_data = p;
    }

    PVector(const PVector & v) = delete;
    PVector& operator=(const PVector& p) = delete;

    static int size() {return DIM;}
    static int dimension() {return DIM;}

    template<class V>
    PVector<F, DIM> & operator += (const V & rhs)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] += rhs[d];
        return *this;
    }

    template<class V>
    PVector<F, DIM> & operator -= (const V & rhs)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] -= rhs[d];
        return *this;
    }

    PVector<F, DIM> & operator *= (const F w)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] *= w;
        return *this;
    }

    PVector<F, DIM> & operator /= (const F w)
    {
        for(auto d = 0; d < DIM; d++)
            m_data[d] /= w;
        return *this;
    }

    F & operator [] (const int i)
    {
        return m_data[i];
    }

    const F & operator [] (const int i) const
    {
        return m_data[i];
    }

private:
    F * m_data;
};

template<typename OS, typename F, int DIM>
OS& operator << (OS & os, const PVector<F, DIM> & p)
{
    if(DIM == 2)
        return os << "Vector_2(" << p[0] << ", " <<p[1] <<')';
    else if(DIM == 3)
        return os << "Vector_3(" << p[0] << ", " <<p[1] <<", " << p[2]<< ')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of PVector_h

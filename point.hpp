#ifndef ASYM_TREE_POINT_HPP
#define ASYM_TREE_POINT_HPP

#include <cstddef>
#include <vector>
#include <cassert>

enum pointCharacterstic
{
    REJECTED_POINT,
    LIVE_POINT
};

template<class _realScalarType>
class point
{
public:
    typedef _realScalarType realScalarType;
    typedef std::vector<realScalarType> realVectorType;
    typedef pointCharacterstic pointCharactersticType;

    point()
    :mCoordinates(1,realScalarType(0)),mDerived(1,realScalarType(0)),mWeight(0),mPointChar(LIVE_POINT)
    {
    }

    point(size_t const numDims,size_t const numDer)
      :mCoordinates(numDims,realScalarType(0))
      ,mDerived(numDer,realScalarType(0))
      ,mWeight()
      ,mPointChar(LIVE_POINT)
    {
    }

    // point(realVectorType const& coordinates,realScalarType const weight)
    // :mCoordinates(coordinates),mWeight(weight),mPointChar(LIVE_POINT),mPointId(0)
    // {
    //     assert(coordinates.size()>0);
    // }

    point(realVectorType const& coordinates,realVectorType const& derived,realScalarType const weight,pointCharactersticType pointChar)
      :mCoordinates(coordinates)
      ,mDerived(derived)
      ,mWeight(weight)
      ,mPointChar(pointChar)
    {
        assert(coordinates.size()>0);
    }

    realScalarType const& weight() const
    {
        return mWeight;
    }

    realScalarType & weight()
    {
        return mWeight;
    }

    size_t size() const
    {
        return mCoordinates.size();
    }

    size_t nder() const
    {
        return mDerived.size();
    }

    realVectorType const & coords() const
    {
        return mCoordinates;
    }

    realVectorType const & derived() const
    {
        return mDerived;
    }

    realScalarType const & operator [] (size_t i) const
    {
        return mCoordinates[i];
    }

    realScalarType & operator [] (size_t i)
    {
        return mCoordinates[i];
    }

    realScalarType const & derived (size_t i) const
    {
        return mDerived[i];
    }

    realScalarType & derived (size_t i)
    {
        return mDerived[i];
    }

    pointCharactersticType const & pointChar() const
    {
        return mPointChar;
    }

    pointCharactersticType & pointChar()
    {
        return mPointChar;
    }

    void set_PointChar(pointCharactersticType new_char)
    {
      mPointChar = new_char;
    }

private:
  realVectorType mCoordinates;
  realVectorType mDerived;
  realScalarType mWeight;
  pointCharactersticType mPointChar;
};

#endif //ASYM_TREE_POINT_HPP

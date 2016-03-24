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

    //static size_t sPointId;

    point()
    :mCoordinates(1,realScalarType(0)),mWeight(0),mPointChar(LIVE_POINT),mPointId(0)
    {
    }

    point(size_t const numDims,realScalarType const weight)
    :mCoordinates(numDims,realScalarType(0)),mWeight(weight),mPointChar(LIVE_POINT),mPointId(0)
    {
    }

    point(realVectorType const& coordinates,realScalarType const weight)
    :mCoordinates(coordinates),mWeight(weight),mPointChar(LIVE_POINT),mPointId(0)
    {
        assert(coordinates.size()>0);
    }

    point(realVectorType const& coordinates,realScalarType const weight,pointCharactersticType pointChar)
    :mCoordinates(coordinates),mWeight(weight),mPointChar(pointChar),mPointId(0)
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

    realVectorType const & coords() const
    {
        return mCoordinates;
    }

    realScalarType const & operator [] (size_t i) const
    {
        return mCoordinates[i];
    }

    realScalarType & operator [] (size_t i)
    {
        return mCoordinates[i];
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

    size_t const & pointId() const
    {
        return mPointId;
    }

    size_t & pointId()
    {
        return mPointId;
    }

private:
    realVectorType mCoordinates;
    realScalarType mWeight;
    pointCharactersticType mPointChar;
    size_t mPointId;
};

// initialise the variable
//template<class _realScalarType > size_t  point<_realScalarType>::sPointId = size_t(0);

#endif //ASYM_TREE_POINT_HPP

#ifndef ASYM_TREE_POINT_HPP
#define ASYM_TREE_POINT_HPP

#include <cstddef>
#include <vector>
#include <cassert>

enum pointCharacterstic
{
    REFERENCE_POINT,
    ACCEPTED_POINT,
    REJECTED_POINT
};

template<class _realScalarType>
class point
{
public:
    typedef _realScalarType realScalarType;
    typedef std::vector<realScalarType> realVectorType;
    typedef pointCharacterstic pointCharactersticType;

    point()
    :mCoordinates(1,realScalarType(0)),mWeight(0),mPointChar(REFERENCE_POINT)
    {

    }

    point(size_t const numDims,realScalarType const weight)
    :mCoordinates(numDims,realScalarType(0)),mWeight(weight),mPointChar(REFERENCE_POINT)
    {

    }

    point(realVectorType const& coordinates,realScalarType const weight)
    :mCoordinates(coordinates),mWeight(weight),mPointChar(REFERENCE_POINT)
    {
        assert(coordinates.size()>0);
    }

    point(realVectorType const& coordinates,realScalarType const weight,pointCharactersticType pointChar)
    :mCoordinates(coordinates),mWeight(weight),mPointChar(pointChar)
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

private:
    realVectorType mCoordinates;
    realScalarType mWeight;
    pointCharactersticType mPointChar;
};

#endif //ASYM_TREE_POINT_HPP

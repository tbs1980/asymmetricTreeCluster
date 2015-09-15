#ifndef ATC_POINT_HPP
#define ATC_POINT_HPP

#include <cstddef>
#include <vector>
#include <cassert>

template<class realScalarType>
class point
{
public:
    typedef std::vector<realScalarType> realVectorType;

    point()
    :mCoordinates(1,realScalarType(0)),mWeight(0)
    {

    }

    point(size_t const numDims,realScalarType const weight)
    :mCoordinates(numDims,realScalarType(0)),mWeight(weight)
    {

    }

    point(realVectorType const& coordinates,realScalarType const weight)
    :mCoordinates(coordinates),mWeight(weight)
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

private:
    realVectorType mCoordinates;
    realScalarType mWeight;
};

#endif//ATC_POINT_HPP

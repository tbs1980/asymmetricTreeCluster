#ifndef ATC_POINT_HPP
#define ATC_POINT_HPP

#include <cstddef>
#include <vector>
#include <cassert>

template<class _realScalarType>
class point
{
public:
    typedef _realScalarType realScalarType;
    typedef std::vector<realScalarType> realVectorType;

    point()
    :mCoordinates(1,realScalarType(0)),mWeight(0),mAccepted(true)
    {

    }

    point(size_t const numDims,realScalarType const weight)
    :mCoordinates(numDims,realScalarType(0)),mWeight(weight),mAccepted(true)
    {

    }

    point(realVectorType const& coordinates,realScalarType const weight)
    :mCoordinates(coordinates),mWeight(weight),mAccepted(true)
    {
        assert(coordinates.size()>0);
    }

    point(realVectorType const& coordinates,realScalarType const weight,bool accepted)
    :mCoordinates(coordinates),mWeight(weight),mAccepted(accepted)
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

    bool accepted() const
    {
        return mAccepted;
    }

private:
    realVectorType mCoordinates;
    realScalarType mWeight;
    bool mAccepted;
};

#endif//ATC_POINT_HPP

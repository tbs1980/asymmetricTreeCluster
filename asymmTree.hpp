#ifndef ASYMM_TREE_HPP
#define ASYMM_TREE_HPP

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>

template<class pointType>
class asymmTree
{
public:
    typedef typename pointType::realScalarType realScalarType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef std::vector<pointType> pointsArrayType;

    static const size_t thresholdForBranching = 100;

    asymmTree()
    :mLeftSubTree(nullptr),mRightSubTree(nullptr),mSplitDimension(0),mMedianIndex(0)
    {

    }

    asymmTree(pointsArrayType  const& points,
        pointType const & boundMin,
        pointType const & boundMax
    )
    :mLeftSubTree(nullptr),mRightSubTree(nullptr),mPoints(points),mSplitDimension(0)
    ,mBoundMin(boundMin),mBoundMax(boundMax)
    {

        typename pointsArrayType::iterator begin = std::begin(mPoints);
        typename pointsArrayType::iterator end = std::end(mPoints);

        /*
        std::cout<<"input point is "<<std::endl;
        for(size_t i=0;i<mPoints.size();++i)
        {
            for(size_t j=0;j<mPoints[i].size();++j)
            {
                std::cout<<mPoints[i][j]<<"\t";
            }
            std::cout<<mPoints[i].weight()<<std::endl;
        }*/

        // for each dimension find the median and the Fisher discriminant
        std::vector<realScalarType> normDiscr(mPoints[0].size());
        for(size_t i=0;i<mPoints[0].size();++i)
        {
            size_t dim = i;

            // sort the points in the current dimension
            std::sort(begin, end,
                [dim]( pointType a, pointType b)
                {
                    return ( a[dim] < b[dim] );
                }
            );

            // median
            /*
            size_t median_ind = mPoints.size()/2;

            while(median_ind != size_t(0) && mPoints[median_ind][dim] == mPoints[size_t(median_ind - 1)][dim] )
            {
                --median_ind; // put all the nodes with equal coord value in the right subtree
            }

            std::cout<<"Median point in dimension "<<dim<<" is "<<mPoints[median_ind][dim]<<std::endl;
            */

            auto rangeSize = std::distance(begin, end);
            auto median = begin + rangeSize/2;

            /*
            while(median != begin && (*median)[dim] == (*(median - 1))[dim] )
            {
                --median;
            }*/

            std::cout<<"Median point in dimension "<<dim<<" is "<<(*median)[dim]<<std::endl;

            // set the new bounds
            pointType boundMinLeft = boundMin;
            pointType boundMaxLeft = (*median);
            pointType boundMinRight = (*median); //TODO is this correct? should this be the next point?
            pointType boundMaxRight = boundMax;

            // find the Euclidian distance between the min and max for the current dimension
            realScalarType distLeft = std::abs( boundMaxLeft[dim] - boundMinLeft[dim] );
            realScalarType distRight =  std::abs( boundMaxRight[dim] - boundMinRight[dim] );

            std::cout<<"left distnace between min and max points in dim "<<dim <<" is "<<distLeft<<"\t"<<boundMaxLeft[dim]<<"\t"<<boundMinLeft[dim]<<std::endl;
            std::cout<<"right distnace between min and max points in dim "<<dim <<" is "<<distRight<<"\t"<<boundMaxRight[dim]<<"\t"<<boundMinRight[dim]<<std::endl;

            // find the lmin and lmax right and left
            auto wMinMaxLeft = std::minmax_element(begin,median,
                [dim]( pointType a, pointType b)
                {
                    return ( a.weight() < b.weight() );
                }
            );
            auto wMinMaxRight = std::minmax_element( (median+1),end,
                [dim]( pointType a, pointType b)
                {
                    return ( a.weight() < b.weight() );
                }
            );

            realScalarType wMinLeftVal = (*wMinMaxLeft.first).weight();
            realScalarType wMaxLeftVal = (*wMinMaxLeft.second).weight();
            realScalarType wMinRightVal = (*wMinMaxRight.first).weight();
            realScalarType wMaxRightVal = (*wMinMaxRight.second).weight();


            std::cout<<"for the left child the min and max liks are "
                <<(*wMinMaxLeft.first).weight()<<"\t"<<(*wMinMaxLeft.second).weight()<<std::endl;
            std::cout<<"for the right child the min and max liks are "
                <<(*wMinMaxRight.first).weight()<<"\t"<<(*wMinMaxRight.second).weight()<<std::endl;

            realScalarType discLeft = std::abs(wMaxLeftVal-wMinLeftVal)/distLeft;
            realScalarType discRight = std::abs(wMaxRightVal-wMinRightVal)/distRight;

            std::cout<<"left discriminant for dim "<<dim<< " is "<< discLeft <<std::endl;
            std::cout<<"right discriminant for dim "<<dim<< " is "<< discRight <<std::endl;
            std::cout<<"diff between discriminants for dim "<<dim<<" is "<<std::abs(discRight - discLeft)<<std::endl;
            normDiscr[i] = std::abs(discRight - discLeft);
            std::cout<<std::endl;
        }

        auto discrMin = std::min_element(std::begin(normDiscr),std::end(normDiscr));

        std::cout<<"the lowest value for the discriminant is from dim "
            <<std::distance(std::begin(normDiscr),discrMin)<<" = "<<(*discrMin)<<std::endl;

        mSplitDimension = std::distance(std::begin(normDiscr),discrMin);
    }

    ~asymmTree()
    {
        if(mLeftSubTree != nullptr)
        {
            delete mLeftSubTree;
        }

        if(mRightSubTree != nullptr)
        {
            delete mRightSubTree;
        }
    }

private:

    asymmTreeType* mLeftSubTree;
    asymmTreeType* mRightSubTree;
    pointsArrayType mPoints;
    size_t mSplitDimension;
    pointType mBoundMin;
    pointType mBoundMax;
    size_t mMedianIndex;
};

#endif //ASYMM_TREE_HPP

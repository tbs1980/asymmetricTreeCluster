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

    asymmTree(pointsArrayType   points,
        pointType  boundMin,
        pointType  boundMax
    )
    :mLeftSubTree(nullptr),mRightSubTree(nullptr),mPoints(points),mSplitDimension(0)
    ,mBoundMin(boundMin),mBoundMax(boundMax),mPointIndices(points.size())
    {
        for(size_t i=0;i<mPointIndices.size();++i)
        {
            mPointIndices[i] = i;
        }

        typename std::vector<size_t>::iterator begin = std::begin(mPointIndices);
        typename std::vector<size_t>::iterator end = std::end(mPointIndices);

        // for each dimension find the median and the Fisher discriminant
        std::vector<realScalarType> normDiscr(mPoints[0].size());
        for(size_t i=0;i<mPoints[0].size();++i)
        {
            size_t dim = i;

            // sort the points in the current dimension
            std::sort(begin, end,
                [this,dim]( size_t a, size_t b)
                {
                    return ( mPoints[a][dim] < mPoints[b][dim] );
                }
            );

            auto rangeSize = std::distance(begin, end);
            auto median = begin + rangeSize/2;

            while(median != begin && mPoints[*(median)][dim] == mPoints[*(median - 1)][dim] )
            {
                --median;
            }

            std::cout<<"Median point in dimension "<<dim<<" is "<<mPoints[*(median)][dim]<<std::endl;

            // set the new bounds
            pointType boundMinLeft = boundMin;
            pointType boundMaxLeft = mPoints[*(median)];
            pointType boundMinRight = mPoints[*(median)]; //TODO is this correct? should this be the next point?
            pointType boundMaxRight = boundMax;

            // find the Euclidian distance between the min and max for the current dimension
            realScalarType distLeft = std::abs( boundMaxLeft[dim] - boundMinLeft[dim] );
            realScalarType distRight =  std::abs( boundMaxRight[dim] - boundMinRight[dim] );

            std::cout<<"left distnace between min and max points in dim "<<dim <<" is "<<distLeft<<"\t"<<boundMaxLeft[dim]<<"\t"<<boundMinLeft[dim]<<std::endl;
            std::cout<<"right distnace between min and max points in dim "<<dim <<" is "<<distRight<<"\t"<<boundMaxRight[dim]<<"\t"<<boundMinRight[dim]<<std::endl;

            // find the lmin and lmax right and left
            auto wMinMaxLeft = std::minmax_element(begin,median,
                [this,dim](  size_t a, size_t b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );
            auto wMinMaxRight = std::minmax_element( (median+1),end,
                [this,dim](  size_t a, size_t b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );


            realScalarType wMinLeftVal = mPoints[*wMinMaxLeft.first].weight();
            realScalarType wMaxLeftVal = mPoints[*wMinMaxLeft.second].weight();
            realScalarType wMinRightVal = mPoints[*wMinMaxRight.first].weight();
            realScalarType wMaxRightVal = mPoints[*wMinMaxRight.second].weight();


            std::cout<<"for the left child the min and max liks are "
                <<wMinLeftVal<<"\t"<<wMaxLeftVal<<std::endl;
            std::cout<<"for the right child the min and max liks are "
                <<wMinRightVal<<"\t"<<wMaxRightVal<<std::endl;

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

        if(mPoints.size()>thresholdForBranching)
        {
            std::cout<<"Number of points in this node is greater than the threshold. Branching now"<<std::endl;

            std::cout<<"Sorting the points in dim "<<mSplitDimension<<std::endl;

            // sort the points in the current dimension
            std::sort(begin, end,
                [this]( size_t a, size_t b)
                {
                    return ( mPoints[a][mSplitDimension] < mPoints[b][mSplitDimension] );
                }
            );
            auto rangeSize = std::distance(begin, end);
            auto median = begin + rangeSize/2;
            std::cout<<"median is at "<<std::distance(begin, median)<<std::endl;
            while(median != begin &&
                mPoints[*(median)][mSplitDimension] == mPoints[*(median - 1)][mSplitDimension] )
            {
                --median;
            }
            std::cout<<"after search median is at "<<std::distance(begin, median)<<std::endl;

            std::cout<<"Median point in dimension "<<mSplitDimension<<" is "<<mPoints[*(median)][mSplitDimension]<<std::endl;

            // set the new bounds
            pointType boundMinLeft = boundMin;
            pointType boundMaxLeft = mPoints[*(median)];
            pointType boundMinRight = mPoints[*(median)]; //TODO is this correct? should this be the nextpoint?
            pointType boundMaxRight = boundMax;

            // make points for the left and right tree
            pointsArrayType pointsLeft(std::distance(begin, median));
            pointsArrayType pointsRight(std::distance((median+1),end));

            std::cout<<"size of left subtree = "<<pointsLeft.size()<<std::endl;
            std::cout<<"size of right subtree = "<<pointsRight.size()<<std::endl;

            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = mPoints[mPointIndices[i]];
            }

            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = mPoints[mPointIndices[i+pointsLeft.size()]];
            }

            //asymmTreeType* past = new asymmTreeType(pointsLeft,boundMinLeft,boundMaxLeft);

            mLeftSubTree = new asymmTreeType(pointsLeft,boundMinLeft,boundMaxLeft);
            mRightSubTree = new asymmTreeType(pointsRight,boundMinRight,boundMaxRight);

            //delete past;
        }
        else
        {
            std::cout<<"\n-----NOT branching------\n"<<std::endl;
        }


    }

    ~asymmTree()
    {
        if(mLeftSubTree != nullptr)
        {
            std::cout<<"==== deleting the left tree"<<std::endl;
            delete mLeftSubTree;
        }

        if(mRightSubTree != nullptr)
        {
            std::cout<<"==== deleting the right tree"<<std::endl;
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
    std::vector<size_t> mPointIndices;
};

#endif //ASYMM_TREE_HPP

#ifndef ASYMM_TREE_HPP
#define ASYMM_TREE_HPP

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>

template<class pointType>
class asymmTree
{
public:
    typedef typename pointType::realScalarType realScalarType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef std::vector<pointType> pointsArrayType;

    //static const size_t mThresholdForBranching = 100;

    asymmTree(pointsArrayType const&  points,
        pointType  const&  boundMin,
        pointType  const&  boundMax,
        size_t const thresholdForBranching,
        size_t const treeIndex,
        size_t const level

    )
    :mLeftSubTree(nullptr)
    ,mRightSubTree(nullptr)
    ,mPoints(points)
    ,mSplitDimension(0)
    ,mBoundMin(boundMin)
    ,mBoundMax(boundMax)
    ,mMedianVal(boundMin)
    ,mThresholdForBranching(thresholdForBranching)
    ,mTreeIndex(treeIndex)
    ,mTreeLevel(level)
    ,mWeightMin(0)
    ,mWeightMax(0)
    ,mHasLeftSubTree(false)
    ,mHasRighSubTree(false)
    {
        //std::cout<<"Tree index = "<<treeIndex<<std::endl;
        //std::cout<<"Tree level  = "<<level<<std::endl;


        buildTree();
    }

    ~asymmTree()
    {
        if(mLeftSubTree != nullptr)
        {
            //std::cout<<"==== deleting the left tree"<<std::endl;
            delete mLeftSubTree;
        }

        if(mRightSubTree != nullptr)
        {
            //std::cout<<"==== deleting the right tree"<<std::endl;
            delete mRightSubTree;
        }
    }

    void dumpTree(std::ofstream & outFile)
    {
        outFile<<mTreeIndex<<",";
        for(size_t i=0;i<mBoundMin.size();++i)
        {
            outFile<<mBoundMin[i]<<",";
        }
        for(size_t i=0;i<mBoundMax.size();++i)
        {
            outFile<<mBoundMax[i]<<",";
        }
        outFile<<mWeightMin<<","<<mWeightMax<<std::endl;
        if(mHasLeftSubTree)
        {
            mLeftSubTree->dumpTree(outFile);
        }
        if(mHasRighSubTree)
        {
            mRightSubTree->dumpTree(outFile);
        }
    }

    void addPoints(pointsArrayType const& points)
    {
        /*
        std::cout<<"The input points are"<<std::endl;
        for(size_t i=0;i<points.size();++i)
        {
            for(size_t j=0;j<points[i].size()-1;++j)
            {
                std::cout<<points[i][j]<<"\t";
            }
            std::cout<<points[i][points[i].size()-1]<<std::endl;
        }*/

        std::cout<<"We are adding "<<points.size()<<" more points to the tree"<<std::endl;
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            // 1. then we need to sort points and pass to subtree
            // 2. pass the points to left of the median in split dimension to left
            // 3. pass the points to the right to the right tree
            std::cout<<"this node as sub trees. hence sorting the points "
                <<" in dim "<<mSplitDimension<<std::endl;

            // since we have enough points we can create new tress
            std::vector<size_t> pointIndices(points.size());
            // set the point indices for sorting
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
            }

            typename std::vector<size_t>::iterator begin = std::begin(pointIndices);
            typename std::vector<size_t>::iterator end = std::end(pointIndices);

            // sort the points in the current dimension
            std::sort(begin, end,
                [this,points]( size_t a, size_t b)
                {
                    return ( points[a][mSplitDimension] < points[b][mSplitDimension] );
                }
            );

            // find the nuber of points left to the median
            std::cout<<"median value in this dimension is "<<mMedianVal[mSplitDimension]<<std::endl;
            size_t numPointsLeft = 0;
            for(size_t i=0;i<pointIndices.size();++i)
            {
                //std::cout<<i<<"\t"<<"\t"<<pointIndices[i]
                //    <<"\t"<<points[ pointIndices[i] ][mSplitDimension]<<std::endl;
                if(points[ pointIndices[i] ][mSplitDimension] >= mMedianVal[mSplitDimension])
                {
                    break;
                }
                ++numPointsLeft;
            }
            size_t numPointsRight = points.size()-numPointsLeft;
            std::cout<<"Number of points left and right are "
                <<numPointsLeft<<"\t"<<numPointsRight<<std::endl;

            pointsArrayType pointsLeft(numPointsLeft);
            pointsArrayType pointsRight(numPointsRight);

            //std::cout<<"size of left subtree = "<<pointsLeft.size()<<std::endl;
            //std::cout<<"size of right subtree = "<<pointsRight.size()<<std::endl;
            //std::cout<<"Points on the left "<<std::endl;
            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = points[pointIndices[i]];
                //std::cout<<i<<"\t"<<pointsLeft[i][mSplitDimension]<<std::endl;
            }
            //std::cout<<"Points on the right "<<std::endl;
            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = points[pointIndices[i+pointsLeft.size()]];
                //std::cout<<i<<"\t"<<pointsRight[i][mSplitDimension]<<std::endl;
            }

            if(mHasLeftSubTree)
            {
                mLeftSubTree->addPoints(pointsLeft);
            }
            if(mHasRighSubTree)
            {
                mRightSubTree->addPoints(pointsRight);
            }

        }
        else
        {
            // we are the end of the tree
            // we need to see if branching threshold is reached
            if(mPoints.size() + points.size() > mThresholdForBranching)
            {
                std::cout<<"++++++++We have passed the threshold and hence branching"<<std::endl;
                // we are branching
                // 1. add the current set of points to the points
                // 2. make tree

                pointsArrayType newPoints(mPoints.size() + points.size());

                std::cout<<"Old points size = "<<mPoints.size()<<std::endl;
                std::cout<<"New points size = "<<points.size()<<std::endl;


                for(size_t i=0;i<mPoints.size();++i)
                {
                    newPoints[i] = mPoints[i];
                }

                for(size_t i=0;i<points.size();++i)
                {
                    newPoints[i+mPoints.size()] = points[i];
                }

                // assign the new points back to the provate member
                mPoints = newPoints;

                buildTree();

            }
            else
            {
                std::cout<<"*******We haven't got enoung points yet. Not branching*****"<<std::endl;
                // we are not ready to branch yet
                // jut push back the points
            }
        }
    }

private:

    void buildTree()
    {
        // check if we have enough points for branching
        if(mPoints.size()>mThresholdForBranching)// and mTreeLevel <2)
        {
            //std::cout<<"*******>Creating sub trees now."<<std::endl;

            // since we have enough points we can create new tress
            std::vector<size_t> pointIndices(mPoints.size());
            // set the point indices for sorting
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
            }

            typename std::vector<size_t>::iterator begin = std::begin(pointIndices);
            typename std::vector<size_t>::iterator end = std::end(pointIndices);

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

                // find the median
                auto rangeSize = std::distance(begin, end);
                auto median = begin + rangeSize/2;

                // TODO is this step necessary?
                while(median != begin && mPoints[*(median)][dim] == mPoints[*(median - 1)][dim] )
                {
                    --median;
                }

                // set the new bounds
                pointType const boundMinLeft = mBoundMin;
                pointType const boundMaxLeft = mPoints[*(median)];
                pointType const boundMinRight = mPoints[*(median)]; //TODO is this correct? should this be the next point?
                pointType const boundMaxRight = mBoundMax;

                // find the Euclidian distance between the min and max for the current dimension
                realScalarType const distLeft = std::abs( boundMaxLeft[dim] - boundMinLeft[dim] );
                realScalarType const distRight =  std::abs( boundMaxRight[dim] - boundMinRight[dim] );

                // find the lmin and lmax right and left
                auto wMinMaxLeft = std::minmax_element(begin,median,
                    [this](  size_t a, size_t b)
                    {
                        return ( mPoints[a].weight() < mPoints[b].weight() );
                    }
                );
                auto wMinMaxRight = std::minmax_element( (median+1),end,
                    [this](  size_t a, size_t b)
                    {
                        return ( mPoints[a].weight() < mPoints[b].weight() );
                    }
                );


                realScalarType const wMinLeftVal = mPoints[*wMinMaxLeft.first].weight();
                realScalarType const wMaxLeftVal = mPoints[*wMinMaxLeft.second].weight();
                realScalarType const wMinRightVal = mPoints[*wMinMaxRight.first].weight();
                realScalarType const wMaxRightVal = mPoints[*wMinMaxRight.second].weight();

                // find the left and right discriminats
                realScalarType discLeft = std::abs(wMaxLeftVal-wMinLeftVal)/distLeft;
                realScalarType discRight = std::abs(wMaxRightVal-wMinRightVal)/distRight;

                // find the difference between the two sides
                normDiscr[i] = std::abs(discRight - discLeft);
            }

            // find the dimension with the lowest discriminant
            //auto discrMin = std::min_element(std::begin(normDiscr),std::end(normDiscr));
            auto discrMin = std::max_element(std::begin(normDiscr),std::end(normDiscr));

            //std::cout<<"the lowest value for the discriminant is from dim "
            //    <<std::distance(std::begin(normDiscr),discrMin)<<" = "<<(*discrMin)<<std::endl;

            // set the split dimension as the one with the lowest discriminant
            mSplitDimension = std::distance(std::begin(normDiscr),discrMin);


            // sort the points in the split dimension
            std::sort(begin, end,
                [this]( size_t a, size_t b)
                {
                    return ( mPoints[a][mSplitDimension] < mPoints[b][mSplitDimension] );
                }
            );
            auto rangeSize = std::distance(begin, end);
            auto median = begin + rangeSize/2;

            // TODO is this step really necessary?
            while(median != begin &&
                mPoints[*(median)][mSplitDimension] == mPoints[*(median - 1)][mSplitDimension] )
            {
                --median;
            }

            //std::cout<<"Median point in dimension "<<mSplitDimension<<" is "<<mPoints[*(median)][mSplitDimension]<<std::endl;

            // set the new bounds
            pointType boundMinLeft = mBoundMin;
            pointType boundMaxLeft = mBoundMax;
            pointType boundMinRight = mBoundMin;
            pointType boundMaxRight = mBoundMax;

            // the split dimension will have a new bound coming from median point
            boundMaxLeft[mSplitDimension] = mPoints[*(median)][mSplitDimension];
            boundMinRight[mSplitDimension] = mPoints[*(median)][mSplitDimension];

            // since we have decided the split dimension we can set the node bounds here
            //mBoundMin = boundMin;
            //mBoundMax = boundMax;
            mMedianVal = mPoints[*(median)];
            //std::cout<<"Setting the bounds for the dimension "
            //    <<mSplitDimension<<" as "<<mBoundMin[mSplitDimension]<<"\t"<<mMedianVal[mSplitDimension]
            //    <<"\t"<<mBoundMax[mSplitDimension]<<std::endl;

            // make points for the left and right tree
            pointsArrayType pointsLeft(std::distance(begin, median));
            pointsArrayType pointsRight(std::distance((median+1),end));

            //std::cout<<"size of left subtree = "<<pointsLeft.size()<<std::endl;
            //std::cout<<"size of right subtree = "<<pointsRight.size()<<std::endl;

            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = mPoints[pointIndices[i]];
            }

            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = mPoints[pointIndices[i+pointsLeft.size()]];
            }


            mLeftSubTree = new asymmTreeType(pointsLeft,boundMinLeft,boundMaxLeft,
                mThresholdForBranching,(2*mTreeIndex+1),(mTreeLevel+1));
            mHasLeftSubTree = true;
            mRightSubTree = new asymmTreeType(pointsRight,boundMinRight,boundMaxRight,
                mThresholdForBranching,(2*mTreeIndex+2),(mTreeLevel+1));
            mHasRighSubTree = true;

            //std::cout<<"Clear the root points from level "<<mTreeLevel<<" and index "<<mTreeIndex<<" node"<<std::endl;
            mPoints.clear();

        }
        else
        {
            //std::cout<<"\n-----NOT branching------\n"<<std::endl;

            // then we need to find the min and maximum likelihoods
            std::vector<size_t> pointIndices(mPoints.size());
            // set the point indices for sorting
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
            }

            typename std::vector<size_t>::iterator begin = std::begin(pointIndices);
            typename std::vector<size_t>::iterator end = std::end(pointIndices);

            // find the lmin and lmax right and left
            auto wMinMax = std::minmax_element(begin,end,
                [this](  size_t a, size_t b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );

            mWeightMin = mPoints[*wMinMax.first].weight();
            mWeightMax = mPoints[*wMinMax.second].weight();

            //std::cout<<"Maximum and minimum weights for level "<<mTreeLevel
            //<<" and index "<<mTreeIndex<<" node are "<<mWeightMin<<"\t"<<mWeightMax<<std::endl;
        }

    }


    asymmTreeType* mLeftSubTree;
    asymmTreeType* mRightSubTree;
    pointsArrayType mPoints;
    size_t mSplitDimension;
    pointType mBoundMin;
    pointType mBoundMax;
    pointType mMedianVal;
    size_t mThresholdForBranching;
    size_t mTreeIndex;
    size_t mTreeLevel;
    realScalarType mWeightMin;
    realScalarType mWeightMax;
    bool mHasLeftSubTree;
    bool mHasRighSubTree;

};

#endif //ASYMM_TREE_HPP

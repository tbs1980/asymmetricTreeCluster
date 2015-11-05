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
#include <random>

template<class pointType>
class asymmTree
{
public:
    typedef typename pointType::realScalarType realScalarType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef std::vector<pointType> pointsArrayType;

    asymmTree()
    :mLeftSubTree(nullptr)
    ,mRightSubTree(nullptr)
    ,mPoints(0)
    ,mNumDims(0)
    ,mSplitDimension(0)
    ,mThresholdForBranching(0)
    ,mTreeIndex(0)
    ,mTreeLevel(0)
    ,mWeightMin(0)
    ,mWeightMax(0)
    ,mHasLeftSubTree(false)
    ,mHasRighSubTree(false)
    ,mTreeActive(true)
    {

    }

    asymmTree(pointsArrayType const&  points,
        pointType  const&  boundMin,
        pointType  const&  boundMax,
        size_t const thresholdForBranching,
        size_t const treeIndex,
        size_t const level,
        size_t const splitDimension=size_t(0)
    )
    :mLeftSubTree(nullptr)
    ,mRightSubTree(nullptr)
    ,mPoints(points)
    ,mNumDims(boundMin.size())
    ,mSplitDimension(splitDimension)
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
    ,mTreeActive(true)
    {
        // sanity checks
        assert(thresholdForBranching>0);

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

        // compute the mean and standard deviation
        computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn);

        // finally build tree
        buildTree();
    }

    void setup(pointsArrayType const&  points,
        pointType  const&  boundMin,
        pointType  const&  boundMax,
        size_t const thresholdForBranching,
        size_t const treeIndex,
        size_t const level
    )
    {
        assert(mHasLeftSubTree == false and mHasRighSubTree == false);

        // TODO we need to check if the tree is already built?
        // if tree exists we need to delete it first
        // a delete tree functionality needed
        mLeftSubTree = nullptr;
        mRightSubTree = nullptr;
        mPoints = points;
        mNumDims = boundMin.size();
        mSplitDimension = size_t(0);
        mBoundMin = boundMin;
        mBoundMax = boundMax;
        mMedianVal = boundMin;
        mThresholdForBranching = thresholdForBranching;
        mTreeIndex = treeIndex;
        mTreeLevel = level;
        mWeightMin = realScalarType(0);
        mWeightMax = realScalarType(0);
        mHasLeftSubTree = false;
        mHasRighSubTree = false;
        mTreeActive = true;

        // sanity checks
        assert(thresholdForBranching>0);
        // TODO more checks required here

        // set the point indices for sorting
        std::vector<size_t> pointIndices(mPoints.size());
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

        // compute the mean  and the standard deviation
        computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn);

        // finally build the tree
        buildTree();
    }

    void setup(pointType  const&  boundMin,
        pointType  const&  boundMax,
        size_t const thresholdForBranching,
        size_t const treeIndex,
        size_t const level
    )
    {
        assert(mHasLeftSubTree == false and mHasRighSubTree == false);

        // TODO we need to check if the tree is already built?
        // if tree exists we need to delete it first
        // a delete tree functionality needed
        mLeftSubTree = nullptr;
        mRightSubTree = nullptr;
        mNumDims = boundMin.size();
        mSplitDimension = size_t(0);
        mBoundMin = boundMin;
        mBoundMax = boundMax;
        mMedianVal = boundMin;
        mThresholdForBranching = thresholdForBranching;
        mTreeIndex = treeIndex;
        mTreeLevel = level;
        mWeightMin = realScalarType(0);
        mWeightMax = realScalarType(0);
        mHasLeftSubTree = false;
        mHasRighSubTree = false;
        mTreeActive = true;

        assert(thresholdForBranching>0);
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

    void dumpTree(std::ofstream & outFile)
    {
        // TODO should we write a header?
        // the order is
        // tree-index,split-dimension,points-size,bound-min,bound-max,
        // weight-min,weight-max,weight-mean,weight-std-dvn,

        outFile<<mTreeIndex<<","<<mSplitDimension<<","<<mPoints.size()<<",";
        for(size_t i=0;i<mBoundMin.size();++i)
        {
            outFile<<mBoundMin[i]<<",";
        }
        for(size_t i=0;i<mBoundMax.size();++i)
        {
            outFile<<mBoundMax[i]<<",";
        }
        outFile<<mWeightMin<<","
            <<mWeightMax<<","
            <<mWeightsMean<<","
            <<mWeightsStdDvn<<std::endl;

        // go inside the trees and do the same
        if(mHasLeftSubTree)
        {
            mLeftSubTree->dumpTree(outFile);
        }
        if(mHasRighSubTree)
        {
            mRightSubTree->dumpTree(outFile);
        }
    }

    void addPoints(pointsArrayType const& points,bool const makeTree = true)
    {
        if(points.size()==0)
        {
            std::cout<<"No points provided :)"<<std::endl;
        }
        else if(mHasLeftSubTree or mHasRighSubTree)
        {
            // 1. then we need to sort points and pass to subtree
            // 2. pass the points to left of the median in split dimension to left
            // 3. pass the points to the right to the right tree

            // since we have enough points we can create new tress
            std::vector<size_t> pointIndices(points.size());

            // set the point indices for sorting
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
                assert(points[i].size() == mBoundMin.size());
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
            size_t numPointsLeft = 0;
            for(size_t i=0;i<pointIndices.size();++i)
            {
                if(points[ pointIndices[i] ][mSplitDimension] >= mMedianVal[mSplitDimension])
                {
                    break;
                }
                ++numPointsLeft;
            }
            size_t numPointsRight = points.size()-numPointsLeft;

            // create points array for the left and the right trees
            pointsArrayType pointsLeft(numPointsLeft);
            pointsArrayType pointsRight(numPointsRight);

            assert(pointsLeft.size() + pointsRight.size() == points.size());

            // assign points to left and right
            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = points[pointIndices[i]];
            }
            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = points[pointIndices[i+pointsLeft.size()]];
            }

            // branch if we have trees and points
            if(mHasLeftSubTree and pointsLeft.size()>0)
            {
                mLeftSubTree->addPoints(pointsLeft);
            }
            else if(mHasRighSubTree and pointsRight.size()>0)
            {
                mRightSubTree->addPoints(pointsRight);
            }
            else
            {
                std::cout<<"This means we are trying to add point to a tree that we deleted"<<std::endl;
                std::cout<<"numPointsLeft = "<<numPointsLeft<<std::endl;
                std::cout<<"numPointsRight = "<<numPointsRight<<std::endl;
                std::cout<<"mHasLeftSubTree = "<<mHasLeftSubTree<<std::endl;
                std::cout<<"mHasRighSubTree = "<<mHasRighSubTree<<std::endl;
                std::cout<<"The coordinates of the points in the splitDimension "<<mSplitDimension<< " are"<<std::endl;
                for(size_t i=0;i<points.size();++i)
                {
                    std::cout<<points[i][mSplitDimension]<<std::endl;
                }
                std::cout<<"The median value is "<<mMedianVal[mSplitDimension]<< std::endl;

                abort();
            }

        }
        else
        {
            // we are the end of the tree
            // we need to see if branching threshold is reached
            assert(mThresholdForBranching>0);

            // make sure that we have enough points and we are build a tree
            if(mPoints.size() + points.size() > mThresholdForBranching and
                makeTree == true)
            {
                // we are branching
                // 1. add the current set of points to the points
                // 2. make tree

                // sort points
                pointsArrayType newPoints(mPoints.size() + points.size());

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

                // sort and store the min and max
                // set the point indices for sorting
                std::vector<size_t> pointIndices(mPoints.size());
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

                // compute the mean and standard deviation
                computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn);

                // finally build tree
                buildTree();

            }
            else if(points.size()>0)
            {
                // we are not ready to branch yet
                // jut push back the points

                // create a points array to store all points
                // TODO is push back a better solution?
                pointsArrayType newPoints(mPoints.size() + points.size());

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

                // sort and store the min and max
                // set the point indices for sorting
                std::vector<size_t> pointIndices(mPoints.size());
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

                // compute the mean and standard deviation
                computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn);
           }
        }
    }

    void findNearestNodes(pointType const& point,
        pointType const& dist,
        std::vector<size_t> & inds)
    {
        // Find a set of points in a cube/brick defined by dist
        // do we have left or right sub-trees?
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(point[mSplitDimension] < mMedianVal[mSplitDimension])
            {
                // then we are goin to the left tree
                if(mHasLeftSubTree)
                {
                    mLeftSubTree->findNearestNodes(point,dist,inds);
                }
            }
            else
            {
                // then we are going to the right subtree
                if(mHasRighSubTree)
                {
                    mRightSubTree->findNearestNodes(point,dist,inds);
                }

            }

            // only go to the right one if we pass the median when distance is added
            if(point[mSplitDimension] >= mMedianVal[mSplitDimension] and point[mSplitDimension] - dist[mSplitDimension] < mMedianVal[mSplitDimension])
            {
                // then we are going to the right subtree
                if(mHasLeftSubTree)
                {
                    mLeftSubTree->findNearestNodes(point,dist,inds);
                }
            }
        }
        else
        {
            // we are at the end node
            inds.push_back(mTreeIndex);
        }
    }


    void findNodeAndItsDimensions(pointType const& point,
        pointType & nodeDims,
        size_t & treeIndex,
        bool & nodeFound)
    {
        // find the node a point belongs to and retrieve its diemnsions
        nodeFound = false;
        assert(mBoundMin.size() == nodeDims.size());

        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree and point[mSplitDimension] < mMedianVal[mSplitDimension])
            {
                mLeftSubTree->findNodeAndItsDimensions(point,nodeDims,treeIndex,nodeFound);
            }
            else if(mHasRighSubTree and point[mSplitDimension] >= mMedianVal[mSplitDimension])
            {
                mRightSubTree->findNodeAndItsDimensions(point,nodeDims,treeIndex,nodeFound);
            }
            else
            {
                std::cout<<"It looks like this point"<<point[0]<<"\t"<<point[1]<<" belongs  to a delted tree"<<std::endl;
            }
        }
        else
        {
            for(size_t i=0;i<mBoundMin.size();++i)
            {
                assert(mBoundMax[i] > mBoundMin[i]);
                nodeDims[i] = mBoundMax[i] - mBoundMin[i];
                nodeDims[i] *= 10;
            }
            nodeFound = true;
            treeIndex = mTreeIndex;
        }

    }

    void findNearestNodes(pointType const& point,
        size_t const nodeDist,
        std::vector<size_t> & inds)
    {
        // find the nearest nodeDist nodes of the point

        // find the node the point belongs to
        pointType nodeDims(point.size(),realScalarType(0));
        bool nodeFound(false);
        size_t nodeIndex(0);
        findNodeAndItsDimensions(point,nodeDims,nodeIndex,nodeFound);

        // check if we found the node
        assert(nodeFound);

        // get the indices all the active nodes
        std::vector<size_t> treeInds;
        getTreeIndices(treeInds);

        // sort them in order
        std::sort(treeInds.begin(), treeInds.end());

        // find the current node and its nearest neighbours
        auto iter = std::find(treeInds.begin(), treeInds.end(),nodeIndex);

        size_t numNodesTaken = 0;
        size_t pos = std::distance(treeInds.begin(),iter);

        // start the ones before the current node
        for(long i=pos;i>=0;--i)
        {
            if(pos - i > nodeDist)
            {
                break;
            }
            inds.push_back(treeInds[i]);
        }

        // then the ones after the current node
        for(size_t i=pos+1;i<treeInds.size();++i)
        {
            if(i-pos > nodeDist)
            {
                break;
            }

            inds.push_back(treeInds[i]);
        }

    }

    void deleteNodes(realScalarType const weightStar)
    {
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree)
            {
                mLeftSubTree->deleteNodes(weightStar);
                if (mLeftSubTree->treeIsActive()==false)
                {
                    delete mLeftSubTree;
                    mLeftSubTree = nullptr;
                    mHasLeftSubTree = false;
                }
            }
            if(mHasRighSubTree)
            {
                mRightSubTree->deleteNodes(weightStar);
                if(mRightSubTree->treeIsActive() == false)
                {
                    delete mRightSubTree;
                    mRightSubTree = nullptr;
                    mHasRighSubTree = false;
                }
            }
        }
        //else if(mWeightMax <= weightStar)
        //else if(mWeightsMean + realScalarType(1.)*mWeightsStdDvn <= weightStar)
        else if(mWeightMax + realScalarType(2.)*mWeightsStdDvn < weightStar)
        {

            mPoints.clear();
            mTreeActive = false;
        }
    }

    realScalarType weightMin() const
    {
        return mWeightMin;
    }

    realScalarType weightMax() const
    {
        return mWeightMax;
    }

    bool hasLeftSubTree() const
    {
        return mHasLeftSubTree;
    }

    bool hasRightSubTree() const
    {
        return mHasRighSubTree;
    }

    bool treeIsActive() const
    {
        return mTreeActive;
    }

    void getTreeIndices(std::vector<size_t> & inds) const
    {

        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree)
            {
                mLeftSubTree->getTreeIndices(inds);
            }

            if(mHasRighSubTree)
            {
                mRightSubTree->getTreeIndices(inds);
            }
        }
        else
        {
            inds.push_back(mTreeIndex);
        }
    }

    void getBounds(pointType & bndMin, pointType & bndMax,size_t const treeIndex)
    {
        if(mTreeIndex == treeIndex)
        {
            bndMin = mBoundMin;
            bndMax = mBoundMax;
            for(size_t i=0;i<bndMin.size();++i)
            {
                assert(bndMin[i] < bndMax[i]);
            }
        }
        else
        {
            if(mHasRighSubTree)
            {
                mRightSubTree->getBounds(bndMin,bndMax,treeIndex);
            }

            if(mHasLeftSubTree)
            {
                mLeftSubTree->getBounds(bndMin,bndMax,treeIndex);
            }
        }
    }

    bool searchTreeIndex(size_t const treeIndex)
    {
        if(mTreeIndex == treeIndex)
        {
            return true;
        }
        else
        {
            if(mHasLeftSubTree)
            {
                return mLeftSubTree->searchTreeIndex(treeIndex);
            }
            else if(mHasRighSubTree)
            {
                return mRightSubTree->searchTreeIndex(treeIndex);
            }
            else
            {
                return false;
            }
        }
    }

    template<class RNGType>
    pointType randomPoint(RNGType & rng)
    {
        // returns a uniformly sampled random point from the active nodes.
        // pick a node uniformly from the available ones
        std::vector<size_t> treeInds;
        getTreeIndices(treeInds);

        assert(treeInds.size()>0);

        // if there is only one node then we just need to pick from this one
        size_t nodeSelected = treeInds[0];

        // if there are more than one, we need a random selection
        std::uniform_int_distribution<size_t> distUniInt(size_t(0), size_t( treeInds.size()-1) );

        if(treeInds.size()>1)
        {
            nodeSelected = treeInds[distUniInt(rng)];
        }

        // find the bounds of the node we want to generate a point from
        pointType boundMin;
        pointType boundMax;
        getBounds(boundMin,boundMax,nodeSelected);

        pointType randPnt(mBoundMin.size(),realScalarType(0));

        std::uniform_real_distribution<> distUniReal;
        for(size_t i=0;i<boundMin.size();++i)
        {
            assert(boundMin[i] < boundMax[i]);
            randPnt[i] = boundMin[i] + (boundMax[i]-boundMin[i])*distUniReal(rng);
        }

        return randPnt;
    }

    template<class RNGType>
    pointType randomPoint(pointType const& point, RNGType & rng)
    {
        std::vector<size_t> treeInds;
        findNearestNodes(point,size_t(20),treeInds);

        assert(treeInds.size()>0);

        // if there is only one node then we just need to pick from this one
        size_t nodeSelected = treeInds[0];


        // if there are more than one, we need a random selection
        std::uniform_int_distribution<size_t> distUniInt(size_t(0), size_t( treeInds.size()-1) );

        if(treeInds.size()>1)
        {
            nodeSelected = treeInds[distUniInt(rng)];
        }

        // find the bounds of the node we want to generate a point from
        pointType boundMin;
        pointType boundMax;
        getBounds(boundMin,boundMax,nodeSelected);

        pointType randPnt(mBoundMin.size(),realScalarType(0));

        std::uniform_real_distribution<> distUniReal;
        for(size_t i=0;i<boundMin.size();++i)
        {
            assert(boundMin[i] < boundMax[i]);
            randPnt[i] = boundMin[i] + (boundMax[i]-boundMin[i])*distUniReal(rng);
        }

        return randPnt;
    }


    size_t treeIndex(void) const
    {
        return mTreeIndex;
    }

private:

    size_t findMaxVarDimension()
    {
        size_t dimWithMaxVar = 0;
        realScalarType varMax(0);
        for(size_t dim=0;dim<mPoints[0].size();++dim)
        {
            realScalarType sum(0);
            realScalarType sum2(0);

            for(size_t i=0;i<mPoints.size();++i)
            {
                sum += mPoints[i][dim];
                sum2 += mPoints[i][dim]*mPoints[i][dim];
            }

            sum /= (realScalarType) mPoints.size();
            sum2 /= (realScalarType) mPoints.size();

            realScalarType var = sum2 - sum*sum;
            assert(var >= realScalarType(0));

            if(dim == 0)
            {
                varMax = var;
                dimWithMaxVar = dim;
            }
            else
            {
                if(var > varMax)
                {
                    varMax = var;
                    dimWithMaxVar = dim;
                }
            }
        }

        return dimWithMaxVar;
    }

    void buildTree()
    {

        // check if we have enough points for branching
        assert(mThresholdForBranching >0);
        if(mPoints.size()>mThresholdForBranching)// and mTreeLevel <2)
        {
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
            /*
            std::vector<realScalarType> normDiscr(mPoints[0].size());
            for(size_t dim=0;dim<mPoints[0].size();++dim)
            {
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
                assert(rangeSize>0);

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

                if((distLeft>0) == false)
                {
                    std::cout<<"We have a problem with the distance"<<std::endl;
                    for(size_t ii =0;ii<mPoints.size();++ii)
                    {
                        std::cout<<ii<<"\t"<<mPoints[ii][mSplitDimension]<<std::endl;
                    }
                }

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
                assert(distLeft>0);
                assert(distRight>0);

                realScalarType discLeft = std::abs(wMaxLeftVal-wMinLeftVal);///distLeft;
                realScalarType discRight = std::abs(wMaxRightVal-wMinRightVal);///distRight;

                // find the difference between the two sides
                normDiscr[dim] = std::abs(discRight - discLeft);
            }

            // find the dimension with the lowest discriminant
            auto discrMax = std::max_element(std::begin(normDiscr),std::end(normDiscr));
            */

            // set the split dimension as the one with the lowest discriminant
            //mSplitDimension = std::distance(std::begin(normDiscr),discrMax);

            //mSplitDimension = findMaxVarDimension();

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

            // set the new bounds
            pointType boundMinLeft = mBoundMin;
            pointType boundMaxLeft = mBoundMax;
            pointType boundMinRight = mBoundMin;
            pointType boundMaxRight = mBoundMax;


            // the split dimension will have a new bound coming from median point
            boundMaxLeft[mSplitDimension] = mPoints[*(median)][mSplitDimension];
            boundMinRight[mSplitDimension] = mPoints[*(median)][mSplitDimension];

            assert(boundMinLeft[mSplitDimension] < boundMaxLeft[mSplitDimension]);
            assert(boundMinRight[mSplitDimension] < boundMaxRight[mSplitDimension]);

            // since we have decided the split dimension we can set the node bounds here
            mMedianVal = mPoints[*(median)];

            // make points for the left and right tree
            pointsArrayType pointsLeft(std::distance(begin, median));
            pointsArrayType pointsRight(std::distance((median),end));

            assert(pointsLeft.size() + pointsRight.size() == mPoints.size());

            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = mPoints[ pointIndices[i] ];
            }

            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = mPoints[ pointIndices[i+pointsLeft.size()] ];
            }

            assert(pointsLeft.size()>0);
            assert(pointsRight.size()>0);

            mLeftSubTree = new asymmTreeType(pointsLeft,boundMinLeft,boundMaxLeft,
                mThresholdForBranching,(2*mTreeIndex+1),(mTreeLevel+1),(mSplitDimension +1) % mNumDims);
            mHasLeftSubTree = true;

            mRightSubTree = new asymmTreeType(pointsRight,boundMinRight,boundMaxRight,
                mThresholdForBranching,(2*mTreeIndex+2),(mTreeLevel+1),(mSplitDimension +1) % mNumDims);
            mHasRighSubTree = true;

            mPoints.clear();
        }
        else
        {
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

            // find the standard deviation weights as well
            computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn);


        }
    }

    void computeMeanStdDvnOfWeights(pointsArrayType const & points,realScalarType & mean, realScalarType & stdDvn)
    {
        realScalarType sum = 0;
        realScalarType sum2 = 0;

        for(size_t i=0;i<points.size();++i)
        {
            realScalarType weightVal = points[i].weight();
            sum += weightVal;
            sum2 += weightVal*weightVal;
        }

        sum /= (realScalarType)points.size();
        sum2 /= (realScalarType)points.size();

        mean  = sum;

        assert(sum2 >= sum*sum);
        stdDvn = std::sqrt(sum2 - sum*sum);
    }


    asymmTreeType* mLeftSubTree;
    asymmTreeType* mRightSubTree;
    pointsArrayType mPoints;
    size_t mNumDims;
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
    bool mTreeActive;
    realScalarType mWeightsStdDvn;
    realScalarType mWeightsMean;
};

#endif //ASYMM_TREE_HPP

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
#include <numeric>

enum nodeCharacterstic
{
    ACCEPTED,
    REJECTED,
    ACCEPTED_AND_REJECTED
};

template<class pointType>
struct nodeInformation
{
    typedef typename pointType::realScalarType realScalarType;

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
    realScalarType mVolume;
    nodeCharacterstic mNodeChar;
    realScalarType mAccRatio;
};

template<class pointType>
class asymmTree
{
public:
    typedef typename pointType::realScalarType realScalarType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef std::vector<pointType> pointsArrayType;
    typedef nodeInformation<pointType> nodeInformationType;

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
    ,mVolume(0)
    ,mNodeChar(ACCEPTED)
    ,mAccRatio(1)
    {

    }

    asymmTree(pointsArrayType const&  points,
        pointType  const&  boundMin,
        pointType  const&  boundMax,
        size_t const thresholdForBranching,
        size_t const treeIndex,
        size_t const level,
        size_t const splitDimension
    )
    :mLeftSubTree(nullptr)
    ,mRightSubTree(nullptr)
    ,mPoints(points)
    ,mNumDims(boundMin.size())
    ,mSplitDimension(splitDimension)
    ,mBoundMin(boundMin)
    ,mBoundMax(boundMax)
    ,mThresholdForBranching(thresholdForBranching)
    ,mTreeIndex(treeIndex)
    ,mTreeLevel(level)
    ,mWeightMin(0)
    ,mWeightMax(0)
    ,mHasLeftSubTree(false)
    ,mHasRighSubTree(false)
    ,mTreeActive(true)
    ,mVolume(0)
    ,mNodeChar(ACCEPTED)
    ,mAccRatio(1)
    {
        computeNodeCharacterstics();
        if(points.size() >0)
        {
            buildTree();
        }
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

    void dumpTree(std::ofstream & outFile) const
    {
        // TODO should we write a header?
        outFile
        << mTreeIndex << ","           // 0 
        << mSplitDimension << ","      // 1
        << mPoints.size() << ","       // 2
        << mWeightMin << ","           // 3
        << mWeightMax << ","           // 4
        << mWeightsMean << ","         // 5
        << mWeightsStdDvn << ","       // 6
        << mVolume << ","              // 7 
        << (int)mNodeChar << ","       // 8
        << mAccRatio << ",";           // 9

        for(size_t i=0;i<mBoundMin.size();++i)
        {
            outFile<<mBoundMin[i]<<",";
        }
        for(size_t i=0;i<mBoundMax.size() -1 ;++i)
        {
            outFile<<mBoundMax[i]<<",";
        }
        outFile<<mBoundMax[mBoundMax.size() -1]<<std::endl;

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

    void dumpTreeForDotRecurse(std::ofstream & outFile) const
    {
        if(mHasLeftSubTree)
        {
            outFile<<mTreeIndex<<" -> "<<mLeftSubTree->treeIndex()<<std::endl;
        }

        if(mHasRighSubTree)
        {
            outFile<<mTreeIndex<<" -> "<<mRightSubTree->treeIndex()<<std::endl;
        }

        if(mHasLeftSubTree)
        {
            mLeftSubTree->dumpTreeForDotRecurse(outFile);
        }

        if(mHasRighSubTree)
        {
            mRightSubTree->dumpTreeForDotRecurse(outFile);
        }
    }

    void dumpTreeForDot(std::ofstream & outFile) const
    {
        outFile<<"digraph graphname {\n";
        dumpTreeForDotRecurse(outFile);
        outFile<<"}";
    }


    void addPoint(pointType const&  point, bool const makeTree = true)
    {
        assert(mThresholdForBranching>0);

        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if( point[mSplitDimension] < mMedianVal[mSplitDimension] )
            {
                if(mHasLeftSubTree)
                {
                    mLeftSubTree->addPoint(point);
                }
                else
                {
                    std::cout<<"This means we are trying to add point to a left-tree that we deleted"<<std::endl;
                    std::cout<<"The current node id = "<<mTreeIndex<<std::endl;
                    abort();
                }
            }
            else
            {
                if(mHasRighSubTree)
                {
                    mRightSubTree->addPoint(point);
                }
                else
                {
                    std::cout<<"This means we are trying to add point to a right-tree that we deleted"<<std::endl;
                    std::cout<<"The current node id = "<<mTreeIndex<<std::endl;
                    abort();
                }
            }
        }
        else
        {
            mPoints.push_back(point);
            computeNodeCharacterstics();

            if(makeTree == true and mPoints.size() + size_t(1) > mThresholdForBranching)
            {
                buildTree();
            }
        }
    }

    void deleteActiveNodeByIndex(size_t const treeIndex)
    {
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasRighSubTree)
            {
                mRightSubTree->deleteActiveNodeByIndex(treeIndex);
                if(mRightSubTree->treeIsActive() == false)
                {
                    delete mRightSubTree;
                    mRightSubTree = nullptr;
                    mHasRighSubTree = false;
                }
            }
            
            if(mHasLeftSubTree)
            {
                mLeftSubTree->deleteActiveNodeByIndex(treeIndex);
                if (mLeftSubTree->treeIsActive()==false)
                {
                    delete mLeftSubTree;
                    mLeftSubTree = nullptr;
                    mHasLeftSubTree = false;
                }
            }
        }
        else if(treeIndex == mTreeIndex)
        {
            mPoints.clear();
            mTreeActive = false;
        }
    }

    void deleteNodes(realScalarType const reductionFactor)
    {
        // step 1 compute the total volume of accepted and accepted-rejected nodes
        std::vector<nodeInformationType> ndInfVect;
        getTreeIndicesAndVolumesAcc(ndInfVect);

        // only proceed if we have more than one node
        if(ndInfVect.size() > 0)
        {
            //std::cout<<"We have "<<ndInfVect.size()<<" acc-nodes "<<std::endl;
            realScalarType accRejVolume = std::accumulate(ndInfVect.begin(), ndInfVect.end(), 
                realScalarType(0),
                [](realScalarType & a, nodeInformationType & b)
                {
                    return a != realScalarType(0) ? a + b.mVolume : b.mVolume ;
                }
                );
            
            //std::cout<<"Total volume of acc nodes "<<accRejVolume<<std::endl;

            // step 2 define alpha
            assert(reductionFactor > realScalarType(0.5) and reductionFactor <= realScalarType(1) );

            // step 3 compute the volume to be reduced if possible
            realScalarType reducedVolume = ( realScalarType(1)/reductionFactor -realScalarType(1) )*accRejVolume;
            
            //std::cout<<"Reduced volume  = "<<reducedVolume<<std::endl;

            // step 4 get the nodes with rejected points
            std::vector<nodeInformationType> ndInfVectRejc;
            getTreeIndicesAndVolumesRejc(ndInfVectRejc);

            if( ndInfVectRejc.size() > 0 )
            {
                //std::cout<<"We have "<<ndInfVectRejc.size()<<" rejec-nodes"<<std::endl;
                // step 5 sort them according to the minimum likelihood of the node
                std::sort(std::begin(ndInfVectRejc),std::end(ndInfVectRejc),
                    [](nodeInformationType const & a, nodeInformationType const & b)
                    {
                        return a.mWeightMin < b.mWeightMin;
                    }
                    );

                // step 6 delete the nodes
                realScalarType volRejc(0);
                for(size_t i=0;i<ndInfVectRejc.size();++i)
                {
                    //std::cout<<i<<"\t"<<ndInfVectRejc[i].mVolume<<"\t"<<volRejc<<std::endl;
                    if(volRejc + ndInfVectRejc[i].mVolume >= reducedVolume)
                    {
                        break;
                    }
                    else
                    {
                        volRejc += ndInfVectRejc[i].mVolume;
                        //std::cout<<"We should delete "<<ndInfVectRejc[i].mTreeIndex<<std::endl;
                        deleteActiveNodeByIndex(ndInfVectRejc[i].mTreeIndex);
                    }
                }
           }
        }
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

    size_t treeIndex() const
    {
        return mTreeIndex;
    }

    nodeInformationType getNodeInformation() const
    {
        nodeInformationType ndInfo;

        ndInfo.mNumDims = mNumDims;
        ndInfo.mBoundMin = mBoundMin;
        ndInfo.mBoundMax = mBoundMax;
        ndInfo.mMedianVal = mMedianVal;
        ndInfo.mThresholdForBranching = mThresholdForBranching;
        ndInfo.mTreeIndex = mTreeIndex;
        ndInfo.mTreeLevel = mTreeLevel;
        ndInfo.mWeightMin = mWeightMin;
        ndInfo.mWeightMax = mWeightMax;
        ndInfo.mHasLeftSubTree = mHasLeftSubTree;
        ndInfo.mHasRighSubTree = mHasRighSubTree;
        ndInfo.mTreeActive = mTreeActive;
        ndInfo.mWeightsStdDvn = mWeightsStdDvn;
        ndInfo.mWeightsMean = mWeightsMean;
        ndInfo.mVolume = mVolume;
        ndInfo.mNodeChar = mNodeChar;
        ndInfo.mAccRatio = mAccRatio;

        return ndInfo;
    }

    void getTreeIndicesAndVolumesAcc(std::vector<nodeInformationType> & nodeInfoVect) const
    {
        // TODO should we have a struct returning all the properties?
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree)
            {
                mLeftSubTree->getTreeIndicesAndVolumesAcc(nodeInfoVect);
            }

            if(mHasRighSubTree)
            {
                mRightSubTree->getTreeIndicesAndVolumesAcc(nodeInfoVect);
            }
        }
        else if(mNodeChar == ACCEPTED) //or mNodeChar == ACCEPTED_AND_REJECTED)
        {
            nodeInfoVect.push_back( getNodeInformation() );
        }
    }

    void getTreeIndicesAndVolumesRejc(std::vector<nodeInformationType> & nodeInfoVect) const
    {
        // TODO should we have a struct returning all the properties?
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree)
            {
                mLeftSubTree->getTreeIndicesAndVolumesRejc(nodeInfoVect);
            }

            if(mHasRighSubTree)
            {
                mRightSubTree->getTreeIndicesAndVolumesRejc(nodeInfoVect);
            }
        }
        else if(mNodeChar == ACCEPTED_AND_REJECTED or mNodeChar == REJECTED)
        {
            nodeInfoVect.push_back( getNodeInformation() );
        }
    }

    void getTreeInformation(std::vector<nodeInformationType> & nodeInfoVect) const
    {
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasLeftSubTree)
            {
                mLeftSubTree->getTreeInformation(nodeInfoVect);
            }

            if(mHasRighSubTree)
            {
                mRightSubTree->getTreeInformation(nodeInfoVect);
            }
        }
        else
        {
            nodeInfoVect.push_back( getNodeInformation() );
        }
    }


    template<class RNGType>
    pointType getRandomPoint(RNGType & rng)
    {
        // step 1 create a sorted list of active nodes by ascending volume
        std::vector<nodeInformationType> ndInfVect;
        getTreeInformation(ndInfVect);

        assert( ndInfVect.size() > size_t(0) );

        std::sort(std::begin(ndInfVect),std::end(ndInfVect),
            [](nodeInformationType const & a, nodeInformationType const & b)
            {
                return a.mVolume < b.mVolume;
            }
            );

        // step 2 store cumulative volumes and reference indices
        std::vector<realScalarType> fcvol(ndInfVect.size());
        fcvol[0] = ndInfVect[0].mVolume;
        for(size_t i=1;i<ndInfVect.size();++i)
        {
            assert(ndInfVect[i].mVolume > realScalarType(0) );
            fcvol[i] = fcvol[i-1] + ndInfVect[i].mVolume;
        }

        // step 3 convert to cumulative fractional volumes
        assert( fcvol[ndInfVect.size()-1] > realScalarType(0) );
        realScalarType icvol = realScalarType(1) / fcvol[ndInfVect.size()-1];
        for(size_t i=0;i<ndInfVect.size();++i)
        {
            fcvol[i] *= icvol;
        }

        // step 4 uniformly select a vloume element 
        // and find the corresponding node 
        std::uniform_real_distribution<> distUniReal;
        realScalarType uniVal = distUniReal(rng);

        auto lowBnd = std::lower_bound(fcvol.begin(),fcvol.end(),uniVal);

        auto idx = std::distance(fcvol.begin(),lowBnd);

        // step 5 generate a random variate from the node bounds
        pointType boundMin = ndInfVect[ idx ].mBoundMin;
        pointType boundMax = ndInfVect[ idx ].mBoundMax;

        pointType randPnt(boundMin.size(),realScalarType(0));
        for(size_t i=0;i<boundMin.size();++i)
        {
            assert(boundMin[i] < boundMax[i]);
            randPnt[i] = boundMin[i] + (boundMax[i]-boundMin[i])*distUniReal(rng);
        }

        return randPnt;
    }

private:

    void computeNodeCharacterstics()
    {
        // sort and store the min and max
        // set the point indices for sorting
        // TODO check if mPoints is already sorted
        assert(mThresholdForBranching > 0);

        // set the volume
        mVolume = realScalarType(1);
        for(size_t i=0;i<mNumDims;++i)
        {
            assert(mBoundMax[i] > mBoundMin[i]);
            mVolume *= mBoundMax[i] - mBoundMin[i];
        }

        if(mPoints.size() > size_t(0))
        {
            std::vector<size_t> pointIndices( mPoints.size() );
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
            }

            // find the lmin and lmax right and left
            auto wMinMax = std::minmax_element(std::begin(pointIndices),std::end(pointIndices),
                [this](  size_t const a, size_t const b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );

            mWeightMin = mPoints[*wMinMax.first].weight();
            mWeightMax = mPoints[*wMinMax.second].weight();

            // compute the mean and standard deviation
            computeMeanStdDvnOfWeights(mPoints,mWeightsMean,mWeightsStdDvn); // TOTO should this be static?


            // flag accept / reject / accept-reject
            size_t numAcc(0);
            size_t numRej(0);
            for(size_t i=0;i<mPoints.size();++i)
            {
                if(mPoints[i].accepted())
                {
                    numAcc += size_t(1);
                }
                else
                {
                    numRej += size_t(1);
                }
            }

            mAccRatio = (realScalarType)numAcc/(realScalarType)mPoints.size();
            
            if(numAcc > size_t(0) and numRej == size_t(0))
            {
                mNodeChar = ACCEPTED;
            }
            else if(numAcc > size_t(0) and numRej > size_t(0))
            {
                mNodeChar = ACCEPTED_AND_REJECTED;
            }
            else if(numAcc == size_t(0) and numRej > size_t(0))
            {
                mNodeChar = REJECTED;
            }
            else if(numAcc == size_t(0) and numRej == size_t(0))
            {
                assert(mPoints.size() == 0); // in this case we should have problem
            }
            else
            {
                std::cout<<"This shold not hapen. numAcc = "<<numAcc<<", numRej="<<numRej<<std::endl;
                abort();
            }
        }
    }

    size_t findMaxVarDimension() const
    {
        size_t dimWithMaxVar = 0;
        realScalarType varMax(0);
        for(size_t dim= size_t(0);dim<mPoints[0].size();++dim)
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

            if(dim == size_t(0))
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

    size_t findMaxFisherInfoDimension() const
    {
        assert(mPoints.size() > size_t(3)); // otherwise we cant use this

        // since we have enough points we can create new tress
        std::vector<size_t> pointIndices(mPoints.size());

        // set the point indices for sorting
        for(size_t i=0;i<pointIndices.size();++i)
        {
            pointIndices[i] = i;
        }

        auto begin = std::begin(pointIndices);
        auto end = std::end(pointIndices);

        // for each dimension find the median and the Fisher discriminant
        std::vector<realScalarType> discrDiff(mPoints[0].size());

        for(size_t dim=size_t(0);dim<mPoints[0].size();++dim)
        {
            // sort the points in the current dimension
            std::sort(begin, end,
                [this,dim]( size_t const a, size_t const b)
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

            // find the lmin and lmax right and left
            auto wMinMaxLeft = std::minmax_element(begin, median,
                [this](  size_t const a, size_t const b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );
            auto wMinMaxRight = std::minmax_element( (median+1), end,
                [this](  size_t const a, size_t const b)
                {
                    return ( mPoints[a].weight() < mPoints[b].weight() );
                }
            );

            realScalarType const wMinLeftVal = mPoints[*wMinMaxLeft.first].weight();
            realScalarType const wMaxLeftVal = mPoints[*wMinMaxLeft.second].weight();
            realScalarType const wMinRightVal = mPoints[*wMinMaxRight.first].weight();
            realScalarType const wMaxRightVal = mPoints[*wMinMaxRight.second].weight();

            // find the left and right discriminats
            realScalarType discrLeft = wMaxLeftVal - wMinLeftVal;//std::abs( wMaxLeftVal - wMinLeftVal );
            realScalarType discrRight = wMaxRightVal - wMinRightVal;//std::abs( wMaxRightVal - wMinRightVal );

            // find the difference between the two sides
            discrDiff[dim] = discrRight - discrLeft;//std::abs( discrRight - discrLeft );

            std::cout<<"dim = "<<dim<<"\t"<<" discr= "<<discrDiff[dim]<<std::endl;
        }

        // find the dimension with the lowest discriminant
        auto discrMax = std::max_element( std::begin(discrDiff),std::end(discrDiff) );

        // set the split dimension as the one with the lowest discriminant
        return (size_t) std::distance( std::begin(discrDiff), discrMax );
    }

    void buildTree()
    {

        // check if we have enough points for branching
        assert(mThresholdForBranching >0);

        if(mPoints.size() > mThresholdForBranching)
        {
            // since we have enough points we can create new tress
            std::vector<size_t> pointIndices(mPoints.size());

            // set the point indices for sorting
            for(size_t i=0;i<pointIndices.size();++i)
            {
                pointIndices[i] = i;
            }

            mSplitDimension = findMaxVarDimension();
            //mSplitDimension = findMaxFisherInfoDimension();

            auto begin = std::begin(pointIndices);
            auto end = std::end(pointIndices);

            // sort the point-indeices in the split dimension
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

            assert( boundMinLeft[mSplitDimension] < boundMaxLeft[mSplitDimension] );
            assert( boundMinRight[mSplitDimension] < boundMaxRight[mSplitDimension] );

            // since we have decided the split dimension we can set the node bounds here
            mMedianVal = mPoints[*(median)];

            // make points for the left and right tree
            pointsArrayType pointsLeft( std::distance(begin, median) );
            pointsArrayType pointsRight( std::distance(median, end) );

            assert( pointsLeft.size() + pointsRight.size() == mPoints.size() );

            for(size_t i=0;i<pointsLeft.size();++i)
            {
                pointsLeft[i] = mPoints[ pointIndices[i] ];
            }

            for(size_t i=0;i<pointsRight.size();++i)
            {
                pointsRight[i] = mPoints[ pointIndices[i+pointsLeft.size()] ];
            }

            assert( pointsLeft.size() > size_t(0) );
            assert( pointsRight.size() > size_t(0) );

            // branch
            mLeftSubTree = new asymmTreeType(pointsLeft,boundMinLeft,boundMaxLeft,
                mThresholdForBranching,(2*mTreeIndex+1),(mTreeLevel+1),
                (mSplitDimension +1) % mNumDims);
            mHasLeftSubTree = true;

            mRightSubTree = new asymmTreeType(pointsRight,boundMinRight,boundMaxRight,
                mThresholdForBranching,(2*mTreeIndex+2),(mTreeLevel+1),
                (mSplitDimension +1) % mNumDims);
            mHasRighSubTree = true;

            // clear the points in the currect branch
            mPoints.clear();

            // TODO set the branch to inactive?
        }
        else
        {
            computeNodeCharacterstics();
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
    realScalarType mVolume;
    nodeCharacterstic mNodeChar;
    realScalarType mAccRatio;
};


#endif //ASYMM_TREE_HPP

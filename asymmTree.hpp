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
    REFERENCE_NODE,
    ACCEPTED_NODE,
    REJECTED_NODE
};

/**
 * \class asymmTree
 * \brief A class for building an asymmetric k-d tree given a set of points
 * \tparam pointType point-type
 *
 * This class builds a k-d tree from a given set of points. The points are
 * always stored in the final nodes, ie, the ones without any sub-nodes. The
 * partitioning crieteria can be chosen at build time. The class contains
 * methods for generating uniform random points from active nodes. The nodes
 * falls outside accepted region can be deleted.
 */
template<class pointType>
class asymmTree
{
public:

    /**
     * \typedef typename pointType::realScalarType realScalarType
     * \brief real floating point type
     */
    typedef typename pointType::realScalarType realScalarType;

    /**
     * \typedef asymmTree<pointType> asymmTreeType
     * \brief define the self type
     */
    typedef asymmTree<pointType> asymmTreeType;

    /**
     * \typedef std::vector<pointType> pointsArrayType;
     * \brief a vector of points
     */
    typedef std::vector<pointType> pointsArrayType;

    /**
     * \brief The default constructor. Allocates no memory and nodes,
     * just the root node is constructed.
     */

    typedef typename pointType::pointCharactersticType pointCharactersticType;

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
    ,mNodeChar(REFERENCE_NODE)
    ,mAccRatio(1)
    ,mLiveMinWeight(std::numeric_limits<realScalarType>::max())
    ,mLiveMaxWeight(std::numeric_limits<realScalarType>::lowest())
    {

    }

    /**
     * \brief A constructor that sets up a nodes using the charactestics
     * specified.
     *
     * \param points points to construct the tree
     * \param boundMin the lower bounds of the node
     * \param boundMax the upper bounds of the node
     * \param thresholdForBranching threshold (number of points) at which the node is split
     * \param treeIndex an index for the node
     * \param level level of the node in the tree
     * \param splitDimension the dimension in which the tree is split
     */
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
    ,mNodeChar(REFERENCE_NODE)
    ,mAccRatio(1)
    {
        computeNodeCharacterstics();
        if(points.size() >0)
        {
            buildTree();
        }
    }

    /**
     * \brief the default destructor
     */
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

    /**
     * \brief A function that dumps the tree information to a file
     * @param outFile a file handle for the output file
     */
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
        for(size_t i=0;i<mBoundMax.size() ;++i)
        {
            outFile<<mBoundMax[i]<<",";
        }

        outFile<<mNodeChar<<std::endl;

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

    /**
     * \brief a function that outputs the tree in a dot diagram format
     * @param outFile a file handle for the output file
     */
    void dumpTreeForDot(std::ofstream & outFile) const
    {
        outFile<<"digraph graphname {\n";
        dumpTreeForDotRecurse(outFile);
        outFile<<"}";
    }

    /**
     * \brief A method that adds a point into the tree
     * @param point    the point to be added
     * @param makeTree a flag to build tree at the same time
     */
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
            if(mHasLeftSubTree and mHasRighSubTree)
            {
              mLiveMinWeight = std::min(mLeftSubTree->liveMinWeight(),mRightSubTree->liveMinWeight());
              mLiveMaxWeight = std::max(mLeftSubTree->liveMaxWeight(),mRightSubTree->liveMaxWeight());
              mVolAccepted   = mLeftSubTree->volAccepted() + mRightSubTree->volAccepted();
              mVolRejected   = mLeftSubTree->volRejected() + mRightSubTree->volRejected();
            }
            else if(mHasLeftSubTree)
            {
              mLiveMinWeight = mLeftSubTree->liveMinWeight();
              mLiveMaxWeight = mLeftSubTree->liveMaxWeight();
              mVolAccepted   = mLeftSubTree->volAccepted();
              mVolRejected   = mLeftSubTree->volRejected();
            }
            else if(mHasRighSubTree)
            {
              mLiveMinWeight = mRightSubTree->liveMinWeight();
              mLiveMaxWeight = mRightSubTree->liveMaxWeight();
              mVolAccepted   = mRightSubTree->volAccepted();
              mVolRejected   = mRightSubTree->volRejected();
            }
            mNodeChar = REFERENCE_NODE;
        }
        else
        {
            auto lb = std::lower_bound(std::begin(mPoints),std::end(mPoints),point,[]( pointType const a, pointType const b)
            {
                return ( a.weight() < b.weight() );
            });
            mPoints.insert(lb,point);

            computeNodeCharacterstics();

            if(makeTree == true and mPoints.size() + size_t(1) > mThresholdForBranching)
            {
                buildTree();
            }
        }
    }

    /**
     * \para A method for deleting a node specified by the index
     * @param treeIndex The tree index of the node to be deleted
     */
    void deleteActiveNodeByIndex(size_t const treeIndex)
    {
        if(mHasLeftSubTree or mHasRighSubTree)
        {
            if(mHasRighSubTree)
            {
                mRightSubTree->deleteActiveNodeByIndex(treeIndex);
                if( mRightSubTree->hasLeftSubTree() == false
                    and mRightSubTree->hasRightSubTree() == false
                    and mRightSubTree->treeIsActive() == false)
                {
                    delete mRightSubTree;
                    mRightSubTree = nullptr;
                    mHasRighSubTree = false;
                }
            }

            if(mHasLeftSubTree)
            {
                mLeftSubTree->deleteActiveNodeByIndex(treeIndex);
                if ( mLeftSubTree->hasLeftSubTree() == false
                    and mLeftSubTree->hasRightSubTree() == false
                    and mLeftSubTree->treeIsActive() == false)
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


   /**
    * \para A method for deleting a node specified by its address
    * @param treeIndex The address of the node to be deleted
    */
   void deleteActiveNodeByPointer(const asymmTreeType * del_node)
   {
     if(mHasLeftSubTree)
     {
       if(mLeftSubTree == del_node)
       {
         mLeftSubTree->mPoints.clear();
         mLeftSubTree->mTreeActive = false;
       }
       else if(del_node->mBoundMin[mSplitDimension] < mMedianVal[mSplitDimension])
       {
         mLeftSubTree->deleteActiveNodeByPointer(del_node);
       }
       if ( mLeftSubTree->hasLeftSubTree() == false
        and mLeftSubTree->hasRightSubTree() == false
        and mLeftSubTree->treeIsActive() == false)
       {
         delete mLeftSubTree;
         mLeftSubTree = nullptr;
         mHasLeftSubTree = false;
       }
     }

     if(mHasRighSubTree)
     {
       if(mRightSubTree == del_node)
       {
         mRightSubTree->mPoints.clear();
         mRightSubTree->mTreeActive = false;
       }
       else if(del_node->mBoundMax[mSplitDimension] > mMedianVal[mSplitDimension])
       {
         mRightSubTree->deleteActiveNodeByPointer(del_node);
       }
       if ( mRightSubTree->hasLeftSubTree() == false
        and mRightSubTree->hasRightSubTree() == false
        and mRightSubTree->treeIsActive() == false)
       {
         delete mRightSubTree;
         mRightSubTree = nullptr;
         mHasRighSubTree = false;
       }
     }

     if(mHasLeftSubTree and mHasRighSubTree)
     {
       mLiveMinWeight = std::min(mLeftSubTree->liveMinWeight(),mRightSubTree->liveMinWeight());
       mLiveMaxWeight = std::max(mLeftSubTree->liveMaxWeight(),mRightSubTree->liveMaxWeight());
       mVolAccepted   = mLeftSubTree->volAccepted() + mRightSubTree->volAccepted();
       mVolRejected   = mLeftSubTree->volRejected() + mRightSubTree->volRejected();
     }
     else if(mHasLeftSubTree)
     {
       mLiveMinWeight = mLeftSubTree->liveMinWeight();
       mLiveMaxWeight = mLeftSubTree->liveMaxWeight();
       mVolAccepted   = mLeftSubTree->volAccepted();
       mVolRejected   = mLeftSubTree->volRejected();
     }
     else if(mHasRighSubTree)
     {
       mLiveMinWeight = mRightSubTree->liveMinWeight();
       mLiveMaxWeight = mRightSubTree->liveMaxWeight();
       mVolAccepted   = mRightSubTree->volAccepted();
       mVolRejected   = mRightSubTree->volRejected();
     }
   }

    /**
     * \brief A function that returns true if left sub tree exisits
     * @return true if left subtree exists, otherwise false
     */
    bool hasLeftSubTree() const
    {
        return mHasLeftSubTree;
    }

    /**
     * \brief A function that returns true if right sub tree exists
     * @return true if right sub tree exists, otherwise false
     */
    bool hasRightSubTree() const
    {
        return mHasRighSubTree;
    }

    /**
     * \brief A function that returns true if tree is active
     * @return true if tree is active, otherwise false
     */
    bool treeIsActive() const
    {
        return mTreeActive;
    }

    /**
     * \brief A function that returns the tree index of the node
     * @return tree index of the node
     */
    size_t treeIndex() const
    {
        return mTreeIndex;
    }

    pointsArrayType const & getPoints() const
    {
        return mPoints;
    }

    /**
     * \brief Build a list of accepted nodes
     * @param nodeInfoVect a vector containing pointers to accepted nodes
     */
    void getTreeIndicesAndVolumesAcc(std::vector<const asymmTreeType *> & nodeInfoVect) const
    {
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
        else if(mNodeChar == ACCEPTED_NODE)
        {
            nodeInfoVect.push_back( this );
        }
    }

    /**
     * \brief Build a list of rejected nodes
     * @param nodeInfoVect a vector containing pointers to rejected nodes
     */
    void getTreeIndicesAndVolumesRejc(std::vector<const asymmTreeType *> & nodeInfoVect) const
    {
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
        else if(mNodeChar == REJECTED_NODE)
        {
            nodeInfoVect.push_back( this );
        }
    }

    /**
     * \brief Build a list of active nodes
     * @param nodeInfoVect a vector containing pointers to all active nodes
     */
    void getTreeInformation(std::vector<const asymmTreeType *> & nodeInfoVect) const
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
        else { nodeInfoVect.push_back( this ); }
    }

    /**
     * \brief A method that build tree.
     */
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
            mTreeActive = false;

            // Minimum likelihood of live points
            mLiveMinWeight = std::min(mLeftSubTree->liveMinWeight(),mRightSubTree->liveMinWeight());
            mLiveMaxWeight = std::max(mLeftSubTree->liveMaxWeight(),mRightSubTree->liveMaxWeight());
            mVolAccepted   = mLeftSubTree->volAccepted() + mRightSubTree->volAccepted();
            mVolRejected   = mLeftSubTree->volRejected() + mRightSubTree->volRejected();
            mNodeChar      = REFERENCE_NODE;
        }
        else
        {
            computeNodeCharacterstics();
        }
    }

    void searchAndReplacePoint(pointType const & point)
    {
        if( mHasLeftSubTree or mHasRighSubTree)
        {
            // go to the left sub tree and point[mSplitDimension] < mMedianVal[mSplitDimension]
            if( point[mSplitDimension] < mMedianVal[mSplitDimension] )
            {
                if(mHasLeftSubTree)
                {
                    mLeftSubTree->searchAndReplacePoint(point);
                }
            }
            else
            {
                if(mHasRighSubTree)
                {
                    mRightSubTree->searchAndReplacePoint(point);
                }
            }
        }
        else
        {
            for(size_t i=0;i<mPoints.size();++i)
            {
               bool replace=true;
               for(size_t j=0;j<mPoints[0].size();++j)
               {
                   if( std::abs( mPoints[i][j] -  point[j] ) > std::numeric_limits<realScalarType>::epsilon() )
                   {
                       replace = false;
                       break;
                   }
               }

               if(replace)
               {
                   mPoints[i] = point;
               }
            }
        }
    }

    pointType       boundMin()       const {return mBoundMin;}
    pointType       boundMax()       const {return mBoundMax;}
    pointType       median()         const {return mMedianVal;}
    realScalarType  volume()         const {return mVolume;}
    realScalarType  weightMax()      const {return mWeightMax;}
    realScalarType  liveMinWeight()  const {return mLiveMinWeight;}
    realScalarType  liveMaxWeight()  const {return mLiveMaxWeight;}
    realScalarType  splitDimension() const {return mSplitDimension;}
    asymmTreeType * leftSubTree()    const {return mLeftSubTree;}
    asymmTreeType * rightSubTree()   const {return mRightSubTree;}
    realScalarType  volAccepted()    const {return mVolAccepted;}
    realScalarType  volRejected()    const {return mVolRejected;}

    void replace_live()
    {
      if(mHasLeftSubTree && mHasRighSubTree && mLeftSubTree->liveMinWeight() == mRightSubTree->liveMinWeight())
      {
        std::cout << "Cannot identify min likelihood live point. Aborting" << std::endl;
        exit(9);
      }
      if(mHasLeftSubTree && mHasRighSubTree && mLeftSubTree->liveMaxWeight() == mRightSubTree->liveMaxWeight())
        {
          std::cout << "Cannot identify max likelihood live point. Aborting" << std::endl;
          exit(7);
        }
      if(mHasLeftSubTree && !mHasRighSubTree)
      {
        mLeftSubTree->replace_live();
        mLiveMinWeight = mLeftSubTree->liveMinWeight();
        mLiveMaxWeight = mLeftSubTree->liveMaxWeight();
        mVolAccepted   = mLeftSubTree->volAccepted();
        mVolRejected   = mLeftSubTree->volRejected();
      }
      else if(!mHasLeftSubTree && mHasRighSubTree)
      {
        mRightSubTree->replace_live();
        mLiveMinWeight = mRightSubTree->liveMinWeight();
        mLiveMaxWeight = mRightSubTree->liveMaxWeight();
        mVolAccepted   = mRightSubTree->volAccepted();
        mVolRejected   = mRightSubTree->volRejected();
      }
      else if(mHasLeftSubTree && mHasRighSubTree)
      {
        if(mLeftSubTree->liveMinWeight() < mRightSubTree->liveMinWeight())
        {
          mLeftSubTree->replace_live();
          mLiveMinWeight = std::min(mLeftSubTree->liveMinWeight(),mRightSubTree->liveMinWeight());
          mLiveMaxWeight = std::max(mLeftSubTree->liveMaxWeight(),mRightSubTree->liveMaxWeight());
        }
        else if(mLeftSubTree->liveMinWeight() > mRightSubTree->liveMinWeight())
        {
          mRightSubTree->replace_live();
          mLiveMinWeight = std::min(mLeftSubTree->liveMinWeight(),mRightSubTree->liveMinWeight());
          mLiveMaxWeight = std::max(mLeftSubTree->liveMaxWeight(),mRightSubTree->liveMaxWeight());
        }
        mVolAccepted   = mLeftSubTree->volAccepted() + mRightSubTree->volAccepted();
        mVolRejected   = mLeftSubTree->volRejected() + mRightSubTree->volRejected();
      }
      else
      {
        auto lowest_live = std::find_if(std::begin(mPoints),std::end(mPoints),[](pointType const a) {return(a.pointChar()==pointCharactersticType::LIVE_POINT);});
        if( lowest_live == std::end(mPoints) )
        {
          mLiveMinWeight = std::numeric_limits<realScalarType>::max();
          mLiveMaxWeight = std::numeric_limits<realScalarType>::lowest();
          mNodeChar      = REJECTED_NODE;
          printf("Error: Could not find minimum likelihood live point. Exiting\n");
          exit(8);
        }
        else
        {
          lowest_live->set_PointChar(pointCharactersticType::REJECTED_POINT);
          lowest_live++;
          if( lowest_live == std::end(mPoints) )
          {
            mLiveMinWeight = std::numeric_limits<realScalarType>::max();
            mLiveMaxWeight = std::numeric_limits<realScalarType>::lowest();
            mNodeChar      = REJECTED_NODE;
          }
          else
          {
            mLiveMinWeight = lowest_live->weight();
            mLiveMaxWeight = (--std::end(mPoints))->weight();
            mNodeChar      = ACCEPTED_NODE;
          }
        }
        if(mNodeChar == ACCEPTED_NODE)
        {
          mVolAccepted = mVolume;
          mVolRejected = 0.0;
        }
        if(mNodeChar == REJECTED_NODE)
        {
          mVolAccepted = 0.0;
          mVolRejected = mVolume;
        }
      }
    }

private:

    /**
     * \brief A function that outputs the tree in a dot diagram format
     * @param outFile a file handle for the output file
     */
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

            std::sort(std::begin(mPoints),std::end(mPoints),[]( pointType const a, pointType const b)
            {
                return ( a.weight() < b.weight() );
            });

            auto lowest_live = std::find_if(std::begin(mPoints),std::end(mPoints),[](pointType const a) {return(a.pointChar()==pointCharactersticType::LIVE_POINT);});
            if( lowest_live == std::end(mPoints) )
            {
              mLiveMinWeight = std::numeric_limits<realScalarType>::max();
              mLiveMaxWeight = std::numeric_limits<realScalarType>::lowest();
              mNodeChar      = REJECTED_NODE;
            }
            else
            {
              mLiveMinWeight = lowest_live->weight();
              mLiveMaxWeight = (--std::end(mPoints))->weight();
              mNodeChar      = ACCEPTED_NODE;
            }
            if(mNodeChar == ACCEPTED_NODE)
            {
              mVolAccepted = mVolume;
              mVolRejected = 0.0;
            }
            if(mNodeChar == REJECTED_NODE)
            {
              mVolAccepted = 0.0;
              mVolRejected = mVolume;
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
    realScalarType mLiveMinWeight;
    realScalarType mLiveMaxWeight;
    realScalarType mVolAccepted;
    realScalarType mVolRejected;
};


#endif //ASYMM_TREE_HPP

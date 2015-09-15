#ifndef ASYMM_TREE_HPP
#define ASYMM_TREE_HPP

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

template<class pointType>
class asymmTree
{
public:
    typedef asymmTree<pointType> asymmTreeType;
    typedef std::vector<pointType> pointsArrayType;

    asymmTree()
    :mLeftSubTree(nullptr),mRightSubTree(nullptr),mSplitDimension(0)
    {

    }

    asymmTree(pointsArrayType  const& points,size_t const splitDimension)
    :mLeftSubTree(nullptr),mRightSubTree(nullptr),mPoints(points),mSplitDimension(0)
    {
        assert(splitDimension < points.size());

        typename pointsArrayType::iterator begin = std::begin(mPoints);
        typename pointsArrayType::iterator end = std::end(mPoints);

        std::cout<<"input was "<<std::endl;
        for(size_t i=0;i<mPoints.size();++i)
        {
            for(size_t j=0;j<mPoints[i].size()-1;++j)
            {
                std::cout<<mPoints[i][j]<<"\t";
            }
            std::cout<<mPoints[i][mPoints[i].size()-1]<<std::endl;
        }
        std::cout<<std::endl;

        for(size_t i=0;i<mPoints[0].size();++i)
        {
            std::cout<<"sorting points in dimension "<<i<<std::endl;

            size_t dim = i;

            std::sort(begin, end,
                [dim]( pointType a, pointType b)
                {
                    return ( a[dim] < b[dim] );
                }
            );

            std::cout<<"sorted output is "<<std::endl;
            for(size_t j=0;j<mPoints.size();++j)
            {
                for(size_t k=0;k<mPoints[j].size()-1;++k)
                {
                    std::cout<<mPoints[j][k]<<"\t";
                }
                std::cout<<mPoints[j][mPoints[j].size()-1]<<std::endl;
            }
            std::cout<<"Point with "
            std::cout<<std::endl;
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

private:

    asymmTreeType* mLeftSubTree;
    asymmTreeType* mRightSubTree;
    pointsArrayType mPoints;
    size_t mSplitDimension;
};

#endif //ASYMM_TREE_HPP

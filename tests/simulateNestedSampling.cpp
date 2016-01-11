#include <asymmTree.hpp>
#include <point.hpp>
#include <random>
#include <fstream>

//#include "treeSampler.hpp"

template<typename realScalarType>
class GaussLikelihood
{
public:
    GaussLikelihood(size_t numDims)
    :mNumDims(numDims)
    {

    }

    realScalarType logLik(std::vector<realScalarType> const & x)
    {
        assert(x.size()==mNumDims);
        // assume that the mean and sigma in all dimensions are (0,1)

        realScalarType llik(0);
        for(size_t i=0;i<mNumDims;++i)
        {
            llik -= x[i]*x[i];
        }

        return 0.5*llik;
    }
private:
    size_t mNumDims;
};

void simulateNS(void)
{
    typedef point<double> pointType;
    typedef std::vector<pointType> pointsArrayType;
    typedef asymmTree<pointType> asymmTreeType;
    //typedef TreeSampler<pointType> asymmTreeType;
    typedef GaussLikelihood<double> GaussLikelihoodType;

    // define the Gauss dist for computing weights
    const size_t numDims = 2;
    GaussLikelihoodType gauss(numDims);

    // define the bounds
    pointType boundMin(numDims,double(0));
    pointType boundMax(numDims,double(0));

    for(size_t i=0;i<numDims;++i)
    {
        boundMin[i] = -double(2);
        boundMax[i] = double(2);
    }

    //boundMin[1] = -2;
    //boundMax[1] = 10;

    pointsArrayType emptyArry;
    size_t threshold(10);
    size_t treeIndex(0);
    size_t level(0);
    size_t splitDimension(0);
    double reductionFactor(0.9);

    // define the tree
    asymmTreeType ast(emptyArry,boundMin,boundMax,threshold,treeIndex,level,splitDimension);

    // create a set of live points
    size_t numLivePoints = 1000;
    pointsArrayType livePoints(numLivePoints);
    std::vector<size_t> livePointInds(numLivePoints);

    std::mt19937 randGen;
    for(size_t i=0;i<numLivePoints;++i)
    {
        livePointInds[i] = i;
        pointType pt = ast.getRandomPoint(randGen);
        double weight = gauss.logLik(pt.coords());
        pt.weight() = weight;
        pt.pointChar() = REFERENCE_POINT;
        livePoints[i] = pt;
        ast.addPoint(pt,false);
        //std::cout<<std::endl;
    }

    //return ;

    std::sort(livePoints.begin(),livePoints.end(),
        [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });

    std::cout<<"At the begining lmin is "<<livePoints[livePointInds[0]].weight()<<std::endl;

    //return;

    // next loop through the sampling process
    size_t numIter = 1;
    for(size_t j=0;j<numIter;++j)
    {
        size_t tot=0;
        size_t acc=0;
        for(size_t i=0;i<numLivePoints;++i)
        {
            //std::cout<<"\n----------Iteration "<<tot<<" ---------------"<<std::endl;

            pointType fpt = livePoints[livePointInds[0]];

            //  get some random points
            pointType pt = ast.getRandomPoint(randGen);

            double weight = gauss.logLik(pt.coords());
            pt.weight() = weight;

            // replace the live lowest point if necessary
            if(pt.weight() > livePoints[livePointInds[0]].weight())
            {
                pt.pointChar() = ACCEPTED_POINT;
                ast.addPoint(pt,true);

                // delete nodes if necessary
                //ast.deleteNodes(reductionFactor);

                // flag the point as rejected
                ast.searchPoint( livePoints[livePointInds[0]] ).pointChar() = REJECTED_POINT;

                // replace the min live point
                livePoints[livePointInds[0]] = pt;

                // increase the acc rate
                ++acc;
            }
            else
            {
                pt.pointChar() = REJECTED_POINT;
                ast.addPoint(pt,true);
            }

            ++tot;

            std::sort(livePoints.begin(),livePoints.end(),
                [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });
        }

        // at the end of each cycle check if are ready to build tree
        //asymmTreeType tempAst(emptyArry,boundMin,boundMax,threshold,treeIndex,level,splitDimension);
        //asymmTreeType tempAst = ast;
        //tempAst.buildTree();
        //size_t numNodesDeleted = tempAst.deleteNodes(reductionFactor);
        //std::cout<<"nodes deleted = "<<numNodesDeleted<<std::endl;
        
    }

    std::ofstream outFile;
    outFile.open("nsTreeIter.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    return;
    
    //accFile.close();

    outFile.open("nsTreeDotPlot.dot",std::ios::trunc);
    ast.dumpTreeForDot(outFile);
    outFile.close();

}

int main(void)
{
    simulateNS();
    return 0;
}

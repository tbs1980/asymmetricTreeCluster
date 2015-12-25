#include <asymmTree.hpp>
#include <point.hpp>
#include <random>
#include <fstream>

#include "treeSampler.hpp"

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

    size_t threshold = 10;
    size_t treeIndex = 0;
    size_t level = 0;

    // define the tree
    asymmTreeType ast;
    pointsArrayType emptyArry;
    ast = asymmTreeType(emptyArry,boundMin,boundMax,threshold,treeIndex,level);

    //return;

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
        livePoints[i] = pt;
        ast.addPoint(pt,false);
        //std::cout<<std::endl;
    }

    //return ;

    std::sort(livePoints.begin(),livePoints.end(),
        [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });

    std::cout<<"At the begining lmin is "<<livePoints[livePointInds[0]].weight()<<std::endl;


    // createa a file for plotting acceptance rates
    //std::ofstream accFile;
    //accFile.open("acceptance2near.dat",std::ios::trunc);
    //accFile.open("acceptance2all.dat",std::ios::trunc);

    // next loop through the sampling process
    size_t numIter = 3000;
    size_t tot=0;
    size_t acc=0;
    for(size_t i=0;i<numIter;++i)
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
            pt.accepted() = true;
            ast.addPoint(pt);

            // delete nodes if necessary
            ast.deleteNodes();

            // replace the min live point
            livePoints[livePointInds[0]] = pt;

            // increase the acc rate
            ++acc;
        }
        else
        {
            pt.accepted() = false;
            ast.addPoint(pt);
        }

        ++tot;

        std::sort(livePoints.begin(),livePoints.end(),
            [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });

        //std::cout<<"lmin is "<<livePoints[livePointInds[0]].weight()<<"\t"<<(double)acc/(double)tot<<std::endl;


        //accFile<<i<<","<<(double)acc/(double)tot<<","<<livePoints[livePointInds[0]].weight()<<std::endl;
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

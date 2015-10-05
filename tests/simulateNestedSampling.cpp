#include <asymmTree.hpp>
#include <point.hpp>
#include <random>
#include <fstream>

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
    typedef GaussLikelihood<double> GaussLikelihoodType;

    // define the Gauss dist for computing weights
    const size_t numDims = 2;
    GaussLikelihoodType gauss(numDims);

    // define the bounds
    pointType boundMin(numDims,double(0));
    pointType boundMax(numDims,double(0));

    for(size_t i=0;i<numDims;++i)
    {
        boundMin[i] = -double(5);
        boundMax[i] = double(5);
    }

    boundMin[1] = -double(5);
    boundMax[1] = double(15);

    size_t threshold = 4;
    size_t treeIndex = 0;
    size_t level = 0;

    // define the tree
    asymmTreeType ast;
    ast.setup(boundMin,boundMax,threshold,treeIndex,level);

    // create a set of live points
    size_t numLivePoints = 1000;
    pointsArrayType livePoints(numLivePoints);
    std::vector<size_t> livePointInds(numLivePoints);

    std::mt19937 randGen;
    for(size_t i=0;i<numLivePoints;++i)
    {
        livePointInds[i] = i;
        pointType pt = ast.randomPoint(randGen);
        double weight = gauss.logLik(pt.coords());
        pt.weight() = weight;
        livePoints[i] = pt;
    }

    auto begin = std::begin(livePointInds);
    auto end = std::end(livePointInds);

    // sort the live points and delete the nodes
    // which have lmax < l*
    std::sort(begin,end,
        [livePoints](size_t const a,size_t const b)
        {
            return livePoints[a].weight() < livePoints[b].weight();
        }
    );

    //td::cout<<"At the begining lmin is "<<livePoints[livePointInds[0]].weight()<<std::endl;

    // add these points to the tree
    ast.addPoints(livePoints);
    // dump
    std::ofstream outFile;
    outFile.open("nsTree.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    // delte the nodes below the lmin
    // this will be the same as we havn't adde any new points
    ast.deleteNodes(livePoints[livePointInds[0]].weight());
    // dump
    outFile.open("nsTreeInitLv.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    // next loop through the sampling process
    size_t numIter = 1000;
    size_t tot=0;
    size_t acc=0;
    for(size_t i=0;i<numIter;++i)
    {
        //  get some random points
        pointType pt = ast.randomPoint(randGen);
        double weight = gauss.logLik(pt.coords());
        pt.weight() = weight;

        // add this to the tree
        pointsArrayType ptArr(1,pt);
        ast.addPoints(ptArr);

        // replace the live lowest point if necessary
        if(pt.weight() > livePoints[livePointInds[0]].weight())
        {
            // delete nodes if necessary
            //std::cout<<"Deleting the nodes below "<<livePoints[livePointInds[0]].weight()<<std::endl;
            ast.deleteNodes(livePoints[livePointInds[0]].weight());

            // replace the min live point
            livePoints[livePointInds[0]] = pt;

            // increase the acc rate
            ++acc;
        }
        ++tot;

        // re-sort the live points and find min live
        std::sort(begin,end,
            [livePoints](size_t const a,size_t const b)
            {
                return livePoints[a].weight() < livePoints[b].weight();
            }
        );
        //std::cout<<"lmin is "<<livePoints[livePointInds[0]].weight()<<std::endl;

        outFile.open("nsTreeIter.dat",std::ios::trunc);
        ast.dumpTree(outFile);
        outFile.close();

        std::cout<<i<<","<<(double)acc/(double)tot<<","<<livePoints[livePointInds[0]].weight()<<std::endl;
    }

}

int main(void)
{
    simulateNS();
    return 0;
}

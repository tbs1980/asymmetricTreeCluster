#include <asymmTree.hpp>
#include <point.hpp>
#include <random>
#include <fstream>
#include <string>

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
        boundMin[i] = -double(5);
        boundMax[i] = double(5);
    }

    //boundMin[1] = -2;
    //boundMax[1] = 10;

    pointsArrayType emptyArry;
    size_t threshold(10);
    size_t treeIndex(0);
    size_t level(0);
    size_t splitDimension(0);
    double reductionFactor(0.6);

    // define the tree
    asymmTreeType ast(emptyArry,boundMin,boundMax,threshold,treeIndex,level,splitDimension);

    // create a set of live points
    size_t numLivePoints = 1000;
    pointsArrayType livePoints(numLivePoints);
    std::vector<size_t> livePointInds(numLivePoints);

    size_t pointId(100);

    std::mt19937 randGen;
    for(size_t i=0;i<numLivePoints;++i)
    {
        livePointInds[i] = i;
        pointType pt = ast.getRandomPoint(randGen);
        double weight = gauss.logLik(pt.coords());
        pt.weight() = weight;
        pt.pointChar() = ACCEPTED_POINT;
        pt.pointId() = pointId;
        ++pointId;
        livePoints[i] = pt;
        ast.addPoint(pt,false);
        //std::cout<<std::endl;
    }

    //return ;

    std::sort(livePoints.begin(),livePoints.end(),
        [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });

    std::cout<<"At the begining lmin is "<<livePoints[livePointInds[0]].weight()<<std::endl;

    //return;
    //std::ofstream outFilePoints;
    //outFilePoints.open("points.dat",std::ios::trunc);

    // next loop through the sampling process
    size_t numIter = 3;
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
            pt.pointId() = pointId;
            ++pointId;

            double weight = gauss.logLik(pt.coords());
            pt.weight() = weight;

            // replace the live lowest point if necessary
            if(pt.weight() > livePoints[livePointInds[0]].weight())
            {

                pt.pointChar() = ACCEPTED_POINT;
                ast.addPoint(pt,true);


                // flag the point that is going to be replaced as rejected in the tree
                //std::cout<<"Searching for point with id "<<livePoints[ livePointInds[0] ].pointId()<<std::endl;
                livePoints[ livePointInds[0] ].pointChar() = REJECTED_POINT;
                ast.searchAndReplacePoint( livePoints[ livePointInds[0] ] ) ;//.pointChar() = REJECTED_POINT;

                // replace the min live point
                livePoints[livePointInds[0]] = pt;
                //std::cout<<"After assigning the point id is "<<livePoints[ livePointInds[0] ].pointId()<<std::endl;


                // delete nodes if necessary
                ast.deleteNodes(reductionFactor);

                // increase the acc rate
                ++acc;
            }
            else
            {
                pt.pointChar() = REJECTED_POINT;
                ast.addPoint(pt,true);
            }

            //outFilePoints<<pt[0]<<"\t"<<pt[1]<<std::endl;

            ++tot;

            std::sort(livePoints.begin(),livePoints.end(),
                [](pointType const & a,pointType const& b){ return a.weight() < b.weight(); });
        }


        // at the end of each cycle check if are ready to build tree
        /*
        asymmTreeType tempAst(emptyArry,boundMin,boundMax,threshold,treeIndex,level,splitDimension);
        pointsArrayType allPoints = ast.getPoints();
        for(size_t i=0;i<allPoints.size();++i)
        {
            tempAst.addPoint(allPoints[i],false);
        }
        //asymmTreeType tempAst = ast;
        tempAst.buildTree();
        size_t numNodesDeleted = tempAst.deleteNodes(reductionFactor);
        std::cout<<"nodes deleted from the temp tree = "<<numNodesDeleted<<std::endl;


        std::cout<<"deleting nodes from the base tree "<<std::endl;
        if(numNodesDeleted > size_t(0))
        {
            ast.buildTree();
            ast.deleteNodes(reductionFactor);
        }*/

        /*
        if(j==0)
        {
            ast.buildTree();
        }*/

        //sast.deleteNodes(reductionFactor);


        std::ofstream outFile;
        std::string fileNameOut("nsTreeIter_");
        fileNameOut = fileNameOut + std::to_string(j)  + std::string(".dat");
        std::cout<<"fileNameOut  = "<<fileNameOut<<std::endl;
        outFile.open(fileNameOut,std::ios::trunc);
        ast.dumpTree(outFile);
        outFile.close();


    }

    //outFilePoints.close();

    std::cout<<"done looping "<<std::endl;

    return ;

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

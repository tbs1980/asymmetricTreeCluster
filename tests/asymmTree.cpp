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

int testConstructor(void)
{
    typedef point<double> pointType;
    typedef std::vector<pointType> pointsArrayType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef GaussLikelihood<double> GaussLikelihoodType;

    std::default_random_engine generator;

    std::uniform_real_distribution<double> distribution(-5.,5.);

    const size_t numDims = 5;
    const size_t numPoints = 1000;

    // define the bounds
    pointType boundMin(numDims,double(0));
    pointType boundMax(numDims,double(0));

    for(size_t i=0;i<numDims;++i)
    {
        boundMin[i] = -double(5);
        boundMax[i] = double(5);
    }

    // define the Gauss dist for computing weights
    GaussLikelihoodType gauss(numDims);

    // define the points array for the tree
    pointsArrayType pts(numPoints);

    for(size_t i=0;i<numPoints;++i)
    {
        // get the coordinates of the point
        std::vector<double> coords(numDims,0.);
        for(size_t j=0;j<numDims;++j)
        {
            coords[j] = distribution(generator);
        }

        // compute the weight
        double weight = gauss.logLik(coords);

        pointType pv(coords,weight);

        pts[i] = pv;
    }

    size_t threshold = 100;
    size_t treeIndex = 0;
    size_t level = 0;

    asymmTreeType ast(pts,boundMin,boundMax,threshold,treeIndex,level);

    std::ofstream outFile;
    outFile.open("tree.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    return EXIT_SUCCESS;
}

int testAsymmTreeInit(void)
{
    typedef point<double> pointType;
    typedef std::vector<pointType> pointsArrayType;
    typedef asymmTree<pointType> asymmTreeType;
    typedef GaussLikelihood<double> GaussLikelihoodType;

    //asymmTreeType ast;

    std::default_random_engine generator;

    std::uniform_real_distribution<double> distribution0(-5.,5.);
    std::uniform_real_distribution<double> distribution1(-5.,15.);

    const size_t numDims = 2;
    const size_t numPoints = 500;

    // define the bounds
    pointType boundMin(numDims,double(0));
    pointType boundMax(numDims,double(0));

    for(size_t i=0;i<numDims;++i)
    {
        if(i==0)
        {
            boundMin[i] = -double(5);
            boundMax[i] = double(5);
        }
        else
        {
            boundMin[i] = -double(5);
            boundMax[i] = double(15);
        }
    }

    // define the Gauss dist for computing weights
    GaussLikelihoodType gauss(numDims);

    // define the points array for the tree
    pointsArrayType pts(numPoints);

    // assign points and their weights
    std::ofstream outFile;
    outFile.open("points.dat",std::ios::trunc);
    for(size_t i=0;i<numPoints;++i)
    {
        // get the coordinates of the point
        std::vector<double> coords(numDims,0.);
        for(size_t j=0;j<numDims;++j)
        {
            if(j==0)
            {
                coords[j] = distribution0(generator);
            }
            else
            {
                coords[j] = distribution1(generator);
            }

        }

        // compute the weight
        double weight = gauss.logLik(coords);


        for(size_t j=0;j<coords.size();++j)
        {
            outFile<<coords[j]<<",";
        }
        outFile<<weight<<std::endl;

        pointType pv(coords,weight);

        pts[i] = pv;
    }
    outFile.close();


    asymmTreeType ast(pts,boundMin,boundMax,200,0,0);

    outFile.open("tree.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    std::vector<size_t> inds;
    ast.getTreeIndices(inds);

    for(size_t i=0;i<inds.size();++i)
    {
        std::cout<<i<<"\t"<<inds[i]<<std::endl;
    }

    //return EXIT_SUCCESS;

    // delete some nodes
    ast.deleteNodes(-50.0);

    outFile.open("delTree.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();

    std::cout<<"After deletion we get"<<std::endl;
    std::vector<size_t> newInds;
    ast.getTreeIndices(newInds);
    for(size_t i=0;i<newInds.size();++i)
    {
        std::cout<<i<<"\t"<<newInds[i]<<std::endl;
    }

    pointType bMin;
    pointType bMax;
    ast.getBounds(bMin,bMax,size_t(3));

    std::cout<<"Bound min is "<<bMin[0]<<"\t"<<bMin[1]<<std::endl;
    std::cout<<"Bound max is "<<bMax[0]<<"\t"<<bMax[1]<<std::endl;

    pointType randPnt = ast.randomPoint();

    outFile.open("newRandPoint.dat",std::ios::trunc);
    for(size_t i=0;i<randPnt.size()-1;++i)
    {
        outFile<<randPnt[i]<<",";
    }
    outFile<<randPnt[randPnt.size()-1]<<std::endl;
    outFile.close();

    return EXIT_SUCCESS;

    // get some more points to add
    outFile.open("newPoints.dat",std::ios::trunc);
    const size_t numNewPoints = 500; //numPoints/2
    pointsArrayType ptsNew(numNewPoints);
    for(size_t i=0;i<numNewPoints;++i)
    {
        // get the coordinates of the point
        std::vector<double> coords(numDims,0.);
        for(size_t j=0;j<numDims;++j)
        {
            if(j==0)
            {
                coords[j] = distribution0(generator);
                outFile<<coords[j]<<",";
            }
            else
            {
                coords[j] = distribution1(generator);
                outFile<<coords[j]<<",";
            }
        }

        // compute the weight
        double weight = gauss.logLik(coords);

        outFile<<weight<<std::endl;

        pointType pv(coords,weight);

        ptsNew[i] = pv;

    }
    outFile.close();

/*
    ast.addPoints(ptsNew);

    outFile.open("newTree.dat",std::ios::trunc);
    ast.dumpTree(outFile);
    outFile.close();
*/
    return EXIT_SUCCESS;
}
int main(void)
{
    int ret = 0;
    ret += (int)testConstructor();
    //ret += (int)testAsymmTreeInit();
    return ret;
}

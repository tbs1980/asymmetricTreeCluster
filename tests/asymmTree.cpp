#include <asymmTree.hpp>
#include <point.hpp>
#include <random>

int testAsymmTreeInit(void)
{
    typedef point<double> pointType;
    typedef std::vector<pointType> pointsArrayType;
    typedef asymmTree<pointType> asymmTreeType;

    //asymmTreeType ast;

    std::default_random_engine generator;

    std::normal_distribution<double> distribution(0.,1.0);

    const size_t numDims = 2;
    const size_t numPoints = 10;

    pointsArrayType pts(numPoints);

    for(size_t i=0;i<numPoints;++i)
    {
        pointType pv(numDims,double(i));
        for(size_t j=0;j<numDims;++j)
        {
            pv[j] = distribution(generator);
        }

        pts[i] = pv;
    }

    asymmTreeType ast(pts,size_t(0));

    return EXIT_SUCCESS;
}
int main(void)
{
    int ret = 0;
    ret += (int)testAsymmTreeInit();
    return ret;
}

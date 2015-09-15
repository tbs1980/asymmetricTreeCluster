#include <point.hpp>
#include <random>
#include <cstdlib>

int testPointInit(void)
{
    typedef double realScalarType;
    typedef point<realScalarType> pointType;
    typedef std::vector<realScalarType> realVectorType;

    const size_t numDims = 10;

    realVectorType coords(10,realScalarType(0.));

    realScalarType weight = 0.9;

    pointType pt(coords,weight);

    if(pt.size() != numDims)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
int main(void)
{
    int ret = 0;
    ret += (int)testPointInit();
    return ret;
}

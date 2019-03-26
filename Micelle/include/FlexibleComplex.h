#ifndef FLEXIBLECOMPLEX_H
#define FLEXIBLECOMPLEX_H

#include <Definitions.h>
#include <PointInfo.h>
#include <Graph.h>

class FlexibleComplex
{
        ST theComplex;

    public:
        FlexibleComplex();
        static void AlphaTest(Alpha_shape_3 & as, std::string filename);
        static void MakePointDistanceMap(Alpha_shape_3 & as, std::map<Point, PointInfo> & pointDistance);
        static MyGraphType OptimizedSpanningTree(Alpha_shape_3 & as, std::string filename);

    protected:
    private:
};

#endif // FLEXIBLECOMPLEX_H

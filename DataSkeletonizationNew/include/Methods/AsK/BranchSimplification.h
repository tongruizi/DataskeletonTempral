#ifndef BRANCHSIMPLIFICATION_H
#define BRANCHSIMPLIFICATION_H
#include "Graph.h"
#include "BranchDetection.h"

class BranchSimplification
{
    public:
        BranchSimplification();
        static void SimplifyIt(std::list<std::list<Point>> & in,
        std::list<std::list<Point>> & out, double epsilon);
        static void Allocation(MyGraphType & G,
        std::list<std::pair<Point, std::list<Point>>> & outputInfo);
        static void PathToGraph(MyGraphType & G, std::list<std::list<Point>> & in);
        static void PathToGraphProper(MyGraphType & G, std::list<std::list<Point>> & in);

    protected:

    private:
};

#endif // BRANCHSIMPLIFICATION_H

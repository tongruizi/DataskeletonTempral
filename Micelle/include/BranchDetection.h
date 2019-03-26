#ifndef BRANCHDETECTION_H
#define BRANCHDETECTION_H
#include "Graph.h"
#include "Definitions.h"

typedef std::pair<vertex_descriptor, MyGraphType*> VertexGraphPair;
typedef std::priority_queue<VertexGraphPair, std::vector<VertexGraphPair>,
        std::function<bool(VertexGraphPair & a, VertexGraphPair & b)>> VertexHeap;
class BranchDetection
{
public:
    BranchDetection();
    static void SimplifyIt(MyGraphType & G, MyGraphType & K, double minvalue, std::string debug, std::string settings);
    static void PointCollector(vertex_descriptor k, MyGraphType & G, std::list<Point> & pointlist);
    static void SimplifyBranchPoints(MyGraphType & G, std::list<std::pair<Point, std::list<Point>>> & outputInfo);

protected:

private:
};

#endif // BRANCHDETECTION_H

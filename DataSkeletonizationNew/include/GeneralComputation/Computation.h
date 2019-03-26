#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <cstdlib>
#include "Definitions.h"
#include "Graph.h"
#include <mlpack/core.hpp>


class Computation
{
public:
    Computation();
    static double unitRandom()
    {
        return ((double) rand() / (RAND_MAX));
    }
    static double distance(MyGraphType & G, vertex_descriptor v1, vertex_descriptor v2)
    {
        return sqrt(CGAL::squared_distance(G[v1].p,G[v2].p));
    }
    static void treeGraph(MyGraphType &G,std::list<boost::graph_traits<MyGraphType>::edge_descriptor> & edges, MyGraphType & Tree);
    static void ComputeDeluanayTriangulation(MyGraphType & G, std::list<Point> & Vector);
    static void computeMST( std::list<Point> & Vector, MyGraphType & savedtree);
    static void MSTSpecialCompute(std::map<Point, MyGraphType::vertex_descriptor> & vertex_map, MyGraphType & savedtree, std::vector<Point> & Vector);
    static void EpsilonSimplification(MyGraphType & G, double epsilon);
    static double AABBError(MyGraphType & G, std::list<Point> & cloud);
    static double AABBDistance(std::list<std::list<Point>> & paths, std::list<Point> & cloud);
    static double MeanSquareErrorAABB(MyGraphType & G, std::list<Point> & cloud);
    static double ConvertPathsToSegments(std::list<std::list<Point>> & paths, std::list<Segment> & segments);








    //! This method should be moved to convertor instead I think


protected:

private:
};

#endif // COMPUTATION_H

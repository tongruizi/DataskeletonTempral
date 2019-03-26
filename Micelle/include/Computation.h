#ifndef COMPUTATION_H
#define COMPUTATION_H
#include<Definitions.h>
#include<Graph.h>
#include"AbstractComplex.h"
//MyGraphType Computation::treeGraph(MyGraphType &G,std::list<boost::graph_traits<MyGraphType>::edge_descriptor> & edges, MyGraphType & Tree)

class Computation
{
    public:
        Computation();
        static MyGraphType computeMST(std::list<Point> & Vector);
        static double distance(MyGraphType & G, vertex_descriptor v1, vertex_descriptor v2);
        static double unitRandom();
        static void AlphaShapeTest(std::list<Point> & pointcloud);
        static size_t BinomialHash(const ST* st, const ST::Simplex_handle & v);
        static size_t BinomialHash(const Chandler & simplexPair );
        static void FiltrationOfAlphaShapes(std::list<Point> & pointcloud, std::string folder);
        static void treeGraph(MyGraphType &G,std::list<boost::graph_traits<MyGraphType>::edge_descriptor> & edges, MyGraphType & Trees);
        static double AABBDistance(std::list<std::list<Point>> & paths, std::list<Point> & cloud);
        static void ComputeDeluanayTriangulation(MyGraphType & G, std::list<Point> & Vector);
        static void BruteNeighborhoodGraph(MyGraphType & G, std::list<Point> pointlist, double epsilon );
        static void EpsilonSimplification(MyGraphType & G, double epsilon);
        static void MSTSpecialCompute(std::map<Point, MyGraphType::vertex_descriptor> & vertex_map, MyGraphType & savedtree, std::vector<Point> & pointlist);
        static double AverageEdgelength(MyGraphType & G);
        static double AABBError(MyGraphType & G, std::list<Point> cloud);


//std::list<Point> & pointcloud
    protected:



    private:
};

#endif // COMPUTATION_H

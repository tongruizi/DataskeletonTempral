#ifndef DUALTREECOMPUTATION_H
#define DUALTREECOMPUTATION_H

#include <mlpack/methods/emst/dtb.hpp>
#include "Graph.h"
#include "Definitions.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core/tree/cover_tree.hpp>
#include "SegmentDistance.h"


//typedef mlpack::neighbor::NearestNeighborSort WWW;
//typedef mlpack::neighbor::NeighborSearchStat<WWW> NST;
//typedef mlpack::tree::CoverTree<SegmentDistance,NST,arma::mat> OurCoverTree;

template <typename TreeType, typename StatType, typename MatType>
using ReducedCoverTree = mlpack::tree::CoverTree<TreeType,StatType,MatType, mlpack::tree::FirstPointIsRoot>;



class DualTreeComputation
{



public:
    DualTreeComputation();
    // virtual ~DualTreeComputation();
    static double ComputeMST(std::list<Point> & p, MyGraphType & G);
    static void ComputeEpsilonNeighborhoodGraph(std::list<Point> & p, MyGraphType & G, double epsilon);
    static void NearestNeighborForTwoKDTrees(arma::Mat<size_t> & results,arma::mat & referenceset, arma::mat & queryset);
    static void NearestNeighborsForLineSegments(std::vector<Segment> & segments, std::list<Point> & pointcloud,arma::Mat<size_t> & results,arma::mat & distances);
    static double ComputeSmallestValueForConnectedComponent(std::list<Point> & p);

protected:

private:
};

#endif // DUALTREECOMPUTATION_H

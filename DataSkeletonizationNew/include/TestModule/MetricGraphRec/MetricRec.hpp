#ifndef METRICREC_HPP_INCLUDED
#define METRICREC_HPP_INCLUDED

#include "Definitions.h"
#include "Graph.h"
#include <mlpack/methods/range_search/range_search.hpp>
#include <unordered_set>
#include "AbstractAlgorithm.h"
class MetricRec : public AbstractAlgorithm
{
    public:
      MetricRec( double r):
        r(r),AbstractAlgorithm("Metric Reconstruction algorithm") {}

        //! The main method (we will cast list to vector later for convinience
        void Run(std::list<Point> & cloudlist,MyGraphType & out);

        //! The other methods:
        void VectorToMat(arma::mat & original, arma::mat & newOne, std::vector<size_t> & indicesSelected);
        void MedialPointComputator(std::vector<int> & indices, arma::mat & matrix, arma::vec & fProduct);
        void Labeling(std::list<Point> & cloudlist, MyGraphType & out);
        void ProcedureForComputingComponents(int cloudnumber, std::unordered_set<int> & branchPoints,
                                     std::vector<std::vector<int>> & pointsToComponents, std::vector<std::vector<size_t>> &resultingNeighbors);
        void ReconstructGraph(std::vector<std::vector<size_t>> &resultingNeighbors, arma::mat & data,
        std::unordered_set<int> & branchPoints, std::unordered_set<int> & edgePoints, MyGraphType & output);


    protected:


    private:
    double r;
    //! Dont store runtime stuff in class:


    //MyGraphType G;
    //std::unordered_set<int> edgePoint;
    //std::unordered_set<int> branchPoint;

};





#endif // METRICREC_HPP_INCLUDED

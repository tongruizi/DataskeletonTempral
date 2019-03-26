#ifndef MORSELAUNCHER_H
#define MORSELAUNCHER_H

#include "AbstractAlgorithm.h"
#include "MorseReconstructionComputator.h"

class MorseLauncher : public AbstractAlgorithm
{
public:
    MorseLauncher(int proportion, double densityE, double parC):
        proportion(proportion),densityE(densityE),parC(parC),AbstractAlgorithm("MorseAlgorithm") {}

    void Run(std::list<Point> & cloudlist, MyGraphType & out)
    {
        arma::mat data;
        GeneralConvertor::ListToMatTransposed(cloudlist,data);
        MorseReconstructionComputator::ComputeWithProportions(out,data,proportion,densityE,parC);
    }
protected:


private:
    int proportion;
    double densityE;
    double parC;
};

#endif // MORSELAUNCHER_H

//class AmstLauncher : public AbstractAlgorithm
//{
//public:
//    AmstLauncher(double densityE, int ExpDensityc1, int ExpDensityc2, double allocatorE):AbstractAlgorithm("Amst"),
//        densityE(densityE), ExpDensityc1(ExpDensityc1), ExpDensityc2(ExpDensityc2), allocatorE(allocatorE) {}
//    void Run(std::list<Point> & cloudlist, MyGraphType & out) override
//    {
//        std::unordered_map<int, std::vector<int>>  clusters;
//        arma::mat in;
//        arma::mat edgeout;
//        int azz = this -> ExpDensityc1;
//        int azz2 = this-> ExpDensityc2;
//        AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, CleverExpDensity<1,4>,mlpack::tree::KDTree> comp;
//        GeneralConvertor::ListToMatTransposed(cloudlist,in);
//        //ClusterAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2, std::unordered_map<int, std::vector<int>> & clusters)
//        comp.ClusterAMSTComputation(in,edgeout,this->densityE,2.0,this->allocatorE,clusters);
//        GeneralConvertor::ArmaMatToGraph(out,edgeout,in);
//
//    }
//
//
//protected:
//
//private:
//    double densityE;
//    int ExpDensityc1;
//    int ExpDensityc2;
//    double allocatorE;
//};

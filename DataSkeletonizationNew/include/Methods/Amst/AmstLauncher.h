#ifndef AMSTLAUNCHER_H
#define AMSTLAUNCHER_H

#include "AbstractAlgorithm.h"
#include "GeneralConvertor.h"
#include "CleverExpDensity.h"
#include "TestPerformer.h"

class AmstLauncher : public AbstractAlgorithm
{
public:
    AmstLauncher(double densityE, int ExpDensityc1, int ExpDensityc2, double allocatorE):AbstractAlgorithm("Amst"),
        densityE(densityE), ExpDensityc1(ExpDensityc1), ExpDensityc2(ExpDensityc2), allocatorE(allocatorE) {}
    void Run(std::list<Point> & cloudlist, MyGraphType & out) override
    {
        std::unordered_map<int, std::vector<int>>  clusters;
        arma::mat in;
        arma::mat edgeout;
        int azz = this -> ExpDensityc1;
        int azz2 = this-> ExpDensityc2;
        AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, CleverExpDensity<1,4>,mlpack::tree::KDTree> comp;
        GeneralConvertor::ListToMatTransposed(cloudlist,in);
        //ClusterAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2, std::unordered_map<int, std::vector<int>> & clusters)
        comp.ClusterAMSTComputation(in,edgeout,this->densityE,2.0,this->allocatorE,clusters);
        GeneralConvertor::ArmaMatToGraph(out,edgeout,in);

    }


protected:

private:
    double densityE;
    int ExpDensityc1;
    int ExpDensityc2;
    double allocatorE;
};

#endif // AMSTLAUNCHER_H

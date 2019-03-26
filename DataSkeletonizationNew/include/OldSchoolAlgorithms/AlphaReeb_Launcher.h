#ifndef ALPHAREEB_LAUNCHER_H
#define ALPHAREEB_LAUNCHER_H

#include "Graph.h"
#include "Definitions.h"
#include "AlphaReeb_Parameters.h"
#include "AlphaReebComputation.h"
#include "Computation.h"
#include "AbstractAlgorithm.h"
#include "DualTreeComputation.h"

class AlphaReeb_Launcher:public AbstractAlgorithm
{
public:
    AlphaReeb_Launcher(AlphaReeb_Parameters & param,double epsilon):
    param(param),epsilon(epsilon),AbstractAlgorithm("alphaReeb") {}
    void Run(std::list<Point> & cloudlist, MyGraphType & out)
    {
        MyGraphType G;
      //  Computation::ComputeDeluanayTriangulation(G, cloudlist);
     //   Computation::EpsilonSimplification(G, epsilon);
        double minvalueRequired = DualTreeComputation::ComputeSmallestValueForConnectedComponent(cloudlist) + 0.1;
        double realparameter = std::max(this->epsilon,minvalueRequired);
        DualTreeComputation::ComputeEpsilonNeighborhoodGraph(cloudlist,G,realparameter);
        MyGraphType Intermediate;
        AlphaReebComputation::Compute(G, param, out, Intermediate);


    }
protected:

private:
    AlphaReeb_Parameters param;
    double epsilon;
};

#endif // ALPHAREEB_LAUNCHER_H

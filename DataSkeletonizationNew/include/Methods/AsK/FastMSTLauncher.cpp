#include "FastMSTLauncher.h"

FastMSTLauncher::FastMSTLauncher():
AbstractAlgorithm("FastMST")
{
    //ctor
}
void FastMSTLauncher::Run(std::list<Point> & cloudlist, MyGraphType & out)
{
    DualTreeComputation::ComputeMST(cloudlist,out);
}



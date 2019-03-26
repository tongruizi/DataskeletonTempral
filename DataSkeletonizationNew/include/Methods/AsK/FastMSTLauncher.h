#ifndef FASTMSTLAUNCHER_H
#define FASTMSTLAUNCHER_H

#include "DualTreeComputation.h"
#include "AbstractAlgorithm.h"

class FastMSTLauncher : public AbstractAlgorithm
{
    public:
        FastMSTLauncher();
        void Run(std::list<Point> & cloudlist, MyGraphType & out) override;

    protected:

    private:
};

#endif // FASTMSTLAUNCHER_H

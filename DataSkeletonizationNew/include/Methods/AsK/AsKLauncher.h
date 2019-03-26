#ifndef ASKLAUNCHER_H
#define ASKLAUNCHER_H

#include "AbstractAlgorithm.h"

class AsKLauncher : public AbstractAlgorithm
{
    public:
        AsKLauncher(double a, double b, double c,std::string settings,std::string name);
        void Run(std::list<Point> & cloudlist, MyGraphType & out) override;

      //  virtual ~AsKLauncher();

    protected:

    private:
    double branch_detection;
    double branchd;
    double approxError;
    double simplificationError;
    int numberOfRuns;
    std::string settings;
};

#endif // ASKLAUNCHER_H

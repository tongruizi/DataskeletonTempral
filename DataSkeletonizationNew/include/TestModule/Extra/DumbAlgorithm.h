#ifndef DUMBALGORITHM_H
#define DUMBALGORITHM_H

#include "AbstractAlgorithm.h"

class DumbAlgorithm : public AbstractAlgorithm
{
    public:
        DumbAlgorithm():AbstractAlgorithm("Dummy") {}
        void Run(std::list<Point> & cloudlist, MyGraphType & out) override
        {


        }

    protected:

    private:
};

#endif // DUMBALGORITHM_H

#ifndef ALPHAREEBCOMPUTATION_H
#define ALPHAREEBCOMPUTATION_H
#include "AlphaReeb_Parameters.h"
#include "Graph.h"


class AlphaReebComputation
{
    public:
    AlphaReebComputation();
    static void Compute(MyGraphType const& input_graph, AlphaReeb_Parameters const& parameters, MyGraphType& alphaReeb_graph, MyGraphType& InterMediate);

    protected:

    private:
};

#endif // ALPHAREEBCOMPUTATION_H

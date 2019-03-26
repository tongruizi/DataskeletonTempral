#ifndef FUNCTIONEVALUATOR_H
#define FUNCTIONEVALUATOR_H

#include "ExponentialDensity.h"

template<class functiontype>
class FunctionEvaluator
{
    public:
        FunctionEvaluator() {}
        static double Evaluate(double k)
        {
        return functiontype::Evaluate(k);
        }
    protected:

    private:
};


#endif // FUNCTIONEVALUATOR_H

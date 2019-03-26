#include "LinearDensity.h"

LinearDensity::LinearDensity()
{
    //ctor
}

double LinearDensity::Evaluate(double d)
{
    if ((d <= 1)&&(d >= 0))
    {
        return 1-d;
    }
    return 0;
}

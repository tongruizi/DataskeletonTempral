#include "ExponentialDensity.h"


ExponentialDensity::ExponentialDensity()
{
    //ctor
}

double ExponentialDensity::Evaluate(double d)
{
return exp(-d);
}

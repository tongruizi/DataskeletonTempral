#ifndef EXPONENTIALDENSITY_H
#define EXPONENTIALDENSITY_H

#include <Includes.h>

class ExponentialDensity
{
    public:
        ExponentialDensity();
        static double Evaluate(double d);
    protected:

    private:
        int c;
};

#endif // EXPONENTIALDENSITY _H

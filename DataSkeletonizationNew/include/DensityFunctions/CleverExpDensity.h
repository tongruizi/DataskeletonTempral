#ifndef CLEVEREXPDENSITY_H
#define CLEVEREXPDENSITY_H
#include <Includes.h>

template<int p,int q>
class CleverExpDensity
{
    public:
        CleverExpDensity() {}
        static double Evaluate(double d)
        {
        double rr = (double) p / (double) q;
        return exp(-d*rr);
        }


    protected:

    private:
};

#endif // CLEVEREXPDENSITY_H

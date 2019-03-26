#ifndef SIMPLEXCOMPARATOR_H
#define SIMPLEXCOMPARATOR_H
#include "Definitions.h"
#include "SimplexDefinitions.h"

class SimplexComparator
{
    public:
        SimplexComparator();
        bool operator()(const Chandler & one,
        const Chandler & two ) const;
    protected:

    private:
};

#endif // SIMPLEXCOMPARATOR_H

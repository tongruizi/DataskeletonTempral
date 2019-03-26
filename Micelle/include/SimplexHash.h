#ifndef SIMPLEXHASH_H
#define SIMPLEXHASH_H

#include "Definitions.h"
#include "SimplexDefinitions.h"

class SimplexHash
{
    public:
        SimplexHash();
        std::size_t operator()(const Chandler & ob) const;
    protected:

    private:
};

#endif // SIMPLEXHASH_H

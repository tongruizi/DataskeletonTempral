#ifndef HASHFORSIMPLEX_H
#define HASHFORSIMPLEX_H

#include "Definitions.h"
#include "binomial_coeff_table.h"

class HashForSimplex
{
binomial_coeff_table hashtable;
    public:
        HashForSimplex(int n, int m);
        std::size_t operator()(const std::vector<int> & ob) const;
        bool operator()(const std::vector<int> & ob, const std::vector<int> & ob1) const;
        std::size_t giveMeHash(const std::vector<int> & ob) const;


    protected:

    private:
};

#endif // HASHFORSIMPLEX_H

#ifndef BINOMIAL_COEFF_TABLE_H
#define BINOMIAL_COEFF_TABLE_H
#include "Definitions.h"

typedef size_t index_t;

class binomial_coeff_table
{
std::vector<std::vector<index_t>> B;
index_t n_max, k_max;
    public:
        binomial_coeff_table(index_t n, index_t k);
        index_t operator()(index_t n, index_t k) const;
    protected:

    private:
};

#endif // BINOMIAL_COEFF_TABLE_H

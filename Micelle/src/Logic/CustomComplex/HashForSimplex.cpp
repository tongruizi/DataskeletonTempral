#include "HashForSimplex.h"

HashForSimplex::HashForSimplex(int n, int m):
    hashtable(n,m)
{
}

std::size_t HashForSimplex::operator()(const std::vector<int> & ob) const
{
    size_t sum = 0;
    int dim =ob.size();
    for (int i = 0; i < dim; i++)
    {
        sum = sum + hashtable(ob[i]-1,dim);
    }
    return sum;
}

std::size_t HashForSimplex::giveMeHash(const std::vector<int> & ob) const
{
    size_t sum = 0;
    int dim =ob.size();
    for (int i = 0; i < dim; i++)
    {
        sum = sum + hashtable(ob[i]-1,dim);
    }
    return sum;
}

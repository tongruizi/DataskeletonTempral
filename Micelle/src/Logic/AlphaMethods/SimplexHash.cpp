#include "SimplexHash.h"
#include "Definitions.h"
#include "Computation.h"
#include "AbstractComplex.h"

SimplexHash::SimplexHash()
{
    //ctor
}

std::size_t SimplexHash::operator()(const Chandler & ob) const
{
return Computation::BinomialHash(ob);
}

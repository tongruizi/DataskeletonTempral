#ifndef EQUALSCOMPLEX_H
#define EQUALSCOMPLEX_H
#include "Definitions.h"

class EqualsComplex
{
    public:
        EqualsComplex();
        bool operator()(const std::vector<int>& lhs, const std::vector<int>& rhs) const;


    protected:

    private:
};

#endif // EQUALSCOMPLEX_H

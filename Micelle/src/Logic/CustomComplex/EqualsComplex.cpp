#include "EqualsComplex.h"

EqualsComplex::EqualsComplex()
{

}

bool EqualsComplex::operator()(const std::vector<int>& lhs, const std::vector<int>& rhs) const
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }
    for (int i = 0; i < lhs.size(); i++)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }

    }
    return true;
}






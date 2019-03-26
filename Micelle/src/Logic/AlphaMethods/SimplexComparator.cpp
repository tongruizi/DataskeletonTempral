#include "SimplexComparator.h"

SimplexComparator::SimplexComparator()
{
    //ctor
}


bool SimplexComparator::operator()(const Chandler & one,
                                   const Chandler & two ) const
{
    for (int i = 0; i < one.size(); i++)
    {
        if (one[i] != two[i])
        {
            return false;
        }
    }
    return true;
}

#ifndef CELLINFO_H
#define CELLINFO_H

#include "Definitions.h"

struct CellInfo
{
    bool boundary;
    bool irreducable;
    bool found;
    double distance;
    Chandler* prev;
};


#endif // CELLINFO_H

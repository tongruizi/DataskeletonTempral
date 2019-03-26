#ifndef SIMPLEXSORT_H
#define SIMPLEXSORT_H

#include "Definitions.h"
#include "CellInfo.h"
#include "SimplexComparator.h"
#include "SimplexHash.h"

class SimplexSort
{
    std::unordered_map<Chandler, CellInfo,
    SimplexHash, SimplexComparator>* infoMap;

    public:
    SimplexSort(std::unordered_map<Chandler, CellInfo,
    SimplexHash, SimplexComparator>* infoMapp);
    bool operator() (Chandler v, Chandler u);


    protected:

    private:
};

#endif //SIMPLEXSORT_H

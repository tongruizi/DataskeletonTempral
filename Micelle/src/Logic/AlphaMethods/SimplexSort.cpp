#include "SimplexSort.h"


SimplexSort::SimplexSort( std::unordered_map<Chandler, CellInfo,
    SimplexHash, SimplexComparator>* infoMapp):
    infoMap(infoMapp)
{

}

bool SimplexSort::operator() (Chandler v, Chandler u)
{
double length1 = ((*(SimplexSort::infoMap)).find(v)->second).distance;
double length2 = ((*(SimplexSort::infoMap)).find(u)->second).distance;
return length1 > length2;
}



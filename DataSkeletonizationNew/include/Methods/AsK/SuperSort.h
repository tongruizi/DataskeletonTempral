#ifndef SUPERSORT_H
#define SUPERSORT_H

#include "Definitions.h"

struct SuperSort
{
    Point p;
    SuperSort(Point p):p(p)
    {

    }
    inline bool operator() (const Point& struct1, const Point& struct2)
    {
        //return (struct1.key < struct2.key);
        return (CGAL::squared_distance(p,struct1)) < (CGAL::squared_distance(p,struct2));
    }

};

#endif // SUPERSORT_H

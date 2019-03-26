#ifndef AABBMETHODS_H
#define AABBMETHODS_H
#include "Definitions.h"
#include "PointInfo.h"

class AABBMethods
{
Tree tree;
    public:
    AABBMethods(std::list<Triangle> & triangles);
    void ConstructTree(std::list<Triangle> & triangles);
    void ComputeDistances(std::list<Point> & pointlist, std::map<Point,PointInfo> & theMap);

    protected:

    private:
};

#endif // AABBMETHODS_H

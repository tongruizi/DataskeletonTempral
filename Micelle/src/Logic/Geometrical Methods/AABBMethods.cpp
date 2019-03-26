#include "AABBMethods.h"
#include "Definitions.h"

AABBMethods::AABBMethods(std::list<Triangle> & triangles):
tree(triangles.begin(),triangles.end())
{
tree.accelerate_distance_queries();
}

void AABBMethods::ConstructTree(std::list<Triangle> & triangles)
{
//   AABBMethods::tree = Tree(triangles.begin(), triangles.end());
}
void AABBMethods::ComputeDistances(std::list<Point> & pointlist, std::map<Point,PointInfo> & theMap)
{
    for (auto it = pointlist.begin(); it != pointlist.end(); it++)
    {
    double kamkam = sqrt(AABBMethods::tree.squared_distance(*it));
    theMap[(*it)].distance = kamkam;
    }

}


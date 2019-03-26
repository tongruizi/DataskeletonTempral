#ifndef STRAIGHTENINGMETHODS_H
#define STRAIGHTENINGMETHODS_H
#include "Definitions.h"
#include "Graph.h"

typedef std::pair<Point,std::list<Point>> PointPlus;

class StraighteningMethods
{
public:
    StraighteningMethods();
    static bool OptimizeSingle(std::vector<PointPlus> & path,  std::list<Point> & thepath, double e);
    static void Optimize(std::list<std::list<Point>> & out, std::vector<std::vector<PointPlus>> & paths, double e);
    static void Allocate(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G, std::list<Point> & cloud);
    static void GraphToPaths(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G);
    static double ClassicStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e);

protected:

private:
};

#endif // STRAIGHTENINGMETHODS_H

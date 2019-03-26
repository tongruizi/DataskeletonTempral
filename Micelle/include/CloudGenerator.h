#ifndef CLOUDGENERATOR_H
#define CLOUDGENERATOR_H

#include "Definitions.h"
#include "Graph.h"

class CloudGenerator
{
    public:
        CloudGenerator();
        static void generatePoints(int n, MyGraphType & G, double epsilon, std::list<Point> & points);
    protected:

    private:
};

#endif // CLOUDGENERATOR_H

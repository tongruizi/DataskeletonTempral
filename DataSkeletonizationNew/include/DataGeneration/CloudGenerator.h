#ifndef CLOUDGENERATOR_H
#define CLOUDGENERATOR_H

#include "Graph.h"
#include <string>
#include <mlpack/core.hpp>
#include <vector>

class CloudGenerator
{
    public:
        CloudGenerator();
        static void generatePoints(int n, MyGraphType & G, double epsilon, std::list<Point> & points);
        static void generatePointsStandartInArmaMat(int gn, int n, double minangle, double scale, double epsilon, arma::mat & pointz);
        static void generatePointsOnGraphInArmaMat(MyGraphType & G, int pointN, double epsilon, arma::mat & data);
    protected:

    private:
};

#endif // CLOUDGENERATOR_H

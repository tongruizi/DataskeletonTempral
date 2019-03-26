#include "CloudGenerator.h"
#include "Computation.h"
#include "Graph.h"
#include "GraphGeneration.h"

CloudGenerator::CloudGenerator()
{
    //ctor
}

Point affineCombo(Point p, Point q, double t)
{
    double x1 = p.x();
    double x2 = q.x();
    double y1 = p.y();
    double y2 = q.y();
    double z1 = p.z();
    double z2 = q.z();
    double x = t*(x1)+(1-t)*(x2);
    double y = t*(y1)+(1-t)*(y2);
    double z = t*(z1)+(1-t)*(z2);
    return Point(x,y,z);
}

Point perturb(Point p, double epsilon)
{
    double changex = epsilon*(-1 + 2*Computation::unitRandom());
    double changey = epsilon*(-1 + 2*Computation::unitRandom());
    double changez = epsilon*(-1 + 2*Computation::unitRandom());
    double x = p.x() + changex;
    double y = p.y() + changey;
    double z = p.z() + changez;
    return Point(x,y,z);

}

void CloudGenerator::generatePoints(int n, MyGraphType & G, double epsilon, std::list<Point> & listt)
{
    std::vector<std::pair<edge_descriptor, double>> edgez;
    auto edgepair = edges(G);
    double sumdistance = 0;
    for (auto it = edgepair.first; it != edgepair.second; it++)
    {
        double d = G[*it].distance;
        edgez.push_back(std::make_pair((*it),sumdistance));
        sumdistance = sumdistance + d;
    }
    for (int i = 0; i < n; i++)
    {
        double rr = Computation::unitRandom();
        double probability = sumdistance * rr;
        int result = -1;
        for (int j = 0; j < edgez.size(); j++)
        {
            if (edgez[j].second > probability )
            {
                result = j - 1;
                break;
            }
        }
        if (result == -1)
        {
            result = edgez.size() - 1;
        }
        edge_descriptor winner = edgez[result].first;
        double t = Computation::unitRandom();
        Point p = G[source(winner,G)].p;
        Point q = G[target(winner,G)].p;
        Point winwinwin = affineCombo(p,q,t);
        Point win = perturb(winwinwin,epsilon);
        listt.push_back(win);

    }
}

void CloudGenerator::generatePointsOnGraphInArmaMat(MyGraphType & G, int pointN, double epsilon, arma::mat & data)
{
    std::list<Point> points;
    CloudGenerator::generatePoints(pointN,G,epsilon,points);
    int k = points.size();
    data.set_size(k,3);
    int j = 0;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        data(j,0) = it->x();
        data(j,1) = it->y();
        data(j,2) = it->z();
        j++;
    }

}

void CloudGenerator::generatePointsStandartInArmaMat(int gn, int n, double minangle, double scale, double epsilon, arma::mat & data)
{
//! Generation
    MyGraphType G;
    GraphGeneration::RandomGraph1(gn,minangle,scale,G);
//! Conversion
    CloudGenerator::generatePointsOnGraphInArmaMat(G,n,epsilon,data);

}



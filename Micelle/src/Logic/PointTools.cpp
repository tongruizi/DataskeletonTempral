#include "PointTools.h"
#include "Definitions.h"

PointTools ::PointTools ()
{
    //ctor
}

Point PointTools::multiplyPoint(double c, Point & p)
{
    double nx = c * p.x();
    double ny = c * p.y();
    double nz = c * p.z();
    return Point(nx,ny,nz);
}
Point PointTools::sumPoint(Point & p, Point & q)
{
    double nx = p.x() + q.x();
    double ny = p.y() + q.y();
    double nz = p.z() + q.z();
    return Point(nx,ny,nz);
}
Point PointTools::findBarycenter(std::vector<Point> & points)
{
    double kumkum = points.size();
    double coefficent = 1 / kumkum;
    Point p = Point(0,0,0);
    for (int i = 0 ; i < points.size(); i++)
    {
        Point q = PointTools::multiplyPoint(coefficent, points[i]);
        p = PointTools::sumPoint(p, q);
    }
    return p;
}
Point PointTools::findBarycenter(std::list<Point> & points)
{
    double bumbum = points.size();
    double coefficent = 1 / bumbum;
  //  std::cout << "DEBUG: " << points.size() << " Coefficent: " << coefficent << std::endl;
    Point p = Point(0,0,0);
    for (auto it = points.begin(); it != points.end(); it++)
    {
        Point q = PointTools::multiplyPoint(coefficent, *it);

        p = PointTools::sumPoint(p, q);
    }
    return p;
}

Point PointTools::convexCombo(Point & p, Point & q, double t)
{
Point a1 = (multiplyPoint(t,p));
Point a2 = multiplyPoint(1-t,q);
return PointTools::sumPoint(a1,a2);
}

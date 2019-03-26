#ifndef POINTTOOLS_H
#define POINTTOOLS_H

#include "Definitions.h"

class PointTools
{
    public:
        PointTools ();
        static Point multiplyPoint(double c, Point & p);
        static Point sumPoint(Point & p, Point & q);
        static Point findBarycenter(std::vector<Point> & points);
        static Point convexCombo(Point & p, Point & q, double t);
        static Point findBarycenter(std::list<Point> & points);

    protected:

    private:
};

#endif POINTTOOLS_H

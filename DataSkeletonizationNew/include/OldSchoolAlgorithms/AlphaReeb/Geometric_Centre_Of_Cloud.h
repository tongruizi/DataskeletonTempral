#pragma once

#include "Data_Pt.h"

Point Geometric_Centre_Of_Cloud ( std::vector<Point>const& cloud )
{
    size_t cloud_size = cloud.size();
    Point geometric_centre = Point(0,0,0);

    for (auto p : cloud)
    {
        geometric_centre += p;
    }

    geometric_centre = geometric_centre / (double)cloud_size;

    return geometric_centre;
}

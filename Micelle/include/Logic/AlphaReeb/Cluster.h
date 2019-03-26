#pragma once

#include "Definitions.h"
#include "Data_Pt.h"

class Cluster
{
    public:

    int interval;
    Point pt;
    std::vector<Data_Pt> cloud;

    Cluster ( std::vector<Data_Pt>const& c, Point p, int i );

    Cluster();
};

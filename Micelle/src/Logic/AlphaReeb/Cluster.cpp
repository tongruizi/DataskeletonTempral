#include "Cluster.h"

Cluster::Cluster ( std::vector<Data_Pt>const& c, Point p, int i )
{
    cloud = c;
    pt = p;
    interval = i;
}

Cluster::Cluster(){}

#pragma once
#include <map>
#include "Definitions.h"

void Distance_Filter ( std::vector<Point> const& cloud, std::multimap<double, int>& filter_multimap )
{
    size_t cloud_size = cloud.size(), index_1 = 0, index_2 = 0;
    double max_dist = 0;

    for (int counter = 1; counter < cloud_size; ++counter)
    {
        double dist = sqrt(CGAL::squared_distance(cloud[counter],cloud[0]));

        if (dist > max_dist)
        {
            max_dist = dist;
            index_1 = counter;
        }
    }

    max_dist = 0;

    for (int counter = 0; counter < cloud_size; ++counter)
    {
        double dist = sqrt(CGAL::squared_distance( cloud[counter],cloud[index_1]));

        if (dist > max_dist)
        {
            max_dist = dist;
            index_2 = counter;
        }
    }

    for (int counter = 0; counter < cloud_size; ++counter)
    {
        double dist = sqrt(CGAL::squared_distance(cloud[counter],cloud[index_2]));

        filter_multimap.insert( std::pair<double, int>(dist, counter) );
    }
}

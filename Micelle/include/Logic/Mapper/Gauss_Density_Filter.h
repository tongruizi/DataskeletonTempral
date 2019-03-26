#pragma once

#include "Definitions.h"

double Gauss_Density ( std::vector<Point>const& cloud, Point const& data_pt, double sigma )
{
    double Gauss_density = 0;

    for (auto p : cloud)
    {
        Gauss_density += exp( -pow( sqrt(CGAL::squared_distance(data_pt, p)), 2 ) / (double)(2 * pow( sigma, 2 )) );
    }

    return Gauss_density;
}

void Gauss_Density_Filter ( std::vector<Point>const& cloud, double sigma, std::multimap<double, int>& filter_multimap )
{
    for (int i = 0; i < cloud.size(); i++)
    {
        double Gauss_density = Gauss_Density( cloud, cloud[i], sigma );

        filter_multimap.insert( std::pair<double, int>(Gauss_density, i) );
    }
}

#pragma once

#include "Definitions.h"

// Mapper.

void IntervalsMapper ( std::multimap<double, int>& filter_multimap, int num_intervals, double overlap_ratio, std::vector<std::pair<double, double>>& interval_endpts )
{
    auto it_1 = filter_multimap.begin();
    double min = it_1->first;
    auto it_2 = filter_multimap.rbegin();
    double max = it_2->first;
    double range = max - min;

    double length_interval = range / (double)(1 + (num_intervals - 1) * (1 - overlap_ratio));

    double interval_start, interval_end;

    for (int counter = 0; counter < num_intervals; ++counter)
    {
        interval_start = min + counter * length_interval * (1 - overlap_ratio);
        interval_end = interval_start + length_interval;

        interval_endpts.push_back( std::make_pair( interval_start, interval_end ) );
    }
}

void Generate_Subclouds ( std::vector<Point>const& cloud, std::multimap<double, int>& filter_multimap, int num_intervals, double overlap_ratio,
                         std::vector<std::vector<Point>>& subcloud, std::vector<std::unordered_map<int,int>> & vertexIntervals)
{
    std::vector<std::pair<double, double>> interval_endpts;
    interval_endpts.clear();

    IntervalsMapper( filter_multimap, num_intervals, overlap_ratio, interval_endpts );

    // Splitting multimap according to interval endpoints.

    std::multimap<double, int>::iterator it_start, it_end;

    std::vector<std::pair<std::multimap<double, int>::iterator, std::multimap<double, int>::iterator>> pointers;
    pointers.clear();

    for (int counter = 0; counter < num_intervals; ++counter)
    {
        it_start = filter_multimap.lower_bound( interval_endpts[counter].first );
        it_end = filter_multimap.upper_bound( interval_endpts[counter].second );

        pointers.push_back(make_pair(it_start, it_end));
    }

    // Assigning points to subclouds.
    for (int counter = 0; counter < num_intervals; ++counter)
    {
        for (auto it = pointers[counter].first; it != pointers[counter].second; ++it)
        {
            subcloud[counter].push_back( cloud[it->second] );
            int indx = subcloud[counter].size() - 1;
            vertexIntervals[it->second].insert({counter, indx});
        }
    }
}


//  for (int counter = 0; counter < num_intervals; ++counter)
//    {
//        for (auto it = pointers[counter].first; it != pointers[counter].second; ++it)
//        {
//            Point QaQ = G[(*it).second].p;
//            vertex_descriptor newGuy = Graph::add_vertex(Graphs[counter],QaQ);
//        //    Graphs[counter][(*it).second].correspondance = *it;
//            vertexIntervals[it->second].insert({counter,newGuy}); // LOOK at this
//            debugPointCloud[counter].push_back(QaQ);
//            //	subcloud[counter].push_back( cloud[it->second] );
//        }
//    }

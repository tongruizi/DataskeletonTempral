#pragma once

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iostream>
#include <map>
#include "Graph.h"

vertex_descriptor FindTheCore(MyGraphType const& g)
{
    int min = 0;
    int max = boost::num_vertices(g);
    int output = min + (rand() % static_cast<int>(max - min + 1));
    std::vector<vertex_descriptor> parents(boost::num_vertices( g ));
    std::vector<double> distances(boost::num_vertices( g ));
    boost::dijkstra_shortest_paths( g, output,
                                    boost::weight_map( boost::get( boost::edge_weight, g ) )
                                    .distance_map( boost::make_iterator_property_map( distances.begin(), boost::get( boost::vertex_index, g ) ) )
                                    .predecessor_map( boost::make_iterator_property_map( parents.begin(), boost::get( boost::vertex_index, g ) ) ) );

    double maxx_distance = 0;
    vertex_descriptor kam = output;

    for (auto it = boost::vertices(g).first ; it != boost::vertices(g).second; it++)
    {
        double curdistance = distances[*it];
        if (curdistance > maxx_distance)
        {
            maxx_distance = curdistance;
            kam = *it;
        }
    }
    return kam;
}

void Dijkstra ( MyGraphType const& g, std::multimap<double, int>& filter_multimap )
{
    std::vector<vertex_descriptor> parents(boost::num_vertices( g ));
    std::vector<double> distances(boost::num_vertices( g ));
    //std::pair<vertex_iter, vertex_iter> VertexPair = boost::vertices( g );
    vertex_descriptor beginning;
    if (boost::num_vertices(g) < 5)
    {
        beginning = 0;
    }
    else
    {
        beginning = FindTheCore(g);
    }

    boost::dijkstra_shortest_paths( g, beginning,
                                    boost::weight_map( boost::get( boost::edge_weight, g ) )
                                    .distance_map( boost::make_iterator_property_map( distances.begin(), boost::get( boost::vertex_index, g ) ) )
                                    .predecessor_map( boost::make_iterator_property_map( parents.begin(), boost::get( boost::vertex_index, g ) ) ) );
    size_t num_vertices = boost::num_vertices( g );

    for (int counter = 0; counter < num_vertices; ++counter)
    {
        filter_multimap.insert( std::pair<double, int>(distances[counter], counter) );
    }
}

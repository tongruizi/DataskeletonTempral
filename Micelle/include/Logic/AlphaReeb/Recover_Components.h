#pragma once

#include "Data_Pt.h"
#include "Graph.h"

void Recover_Components ( MyGraphType const& g, std::vector<int> const& ComponentNumbering, std::vector<MyGraphType>& components)
{
    auto vpair = boost::vertices(g);
    std::vector<int> correspondance(boost::num_vertices(g));

    for (auto it = vpair.first; it != vpair.second; it++)
    {

        Point pp = g[*it].p;
        vertex_descriptor nv = Graph::add_vertex(components[ComponentNumbering[*it]],pp);
        correspondance[*it] = nv;
    }
    auto epair = boost::edges(g);
    for (auto eit = epair.first; eit != epair.second; eit++)
    {
        vertex_descriptor first = boost::source(*eit, g);
        vertex_descriptor second = boost::target(*eit, g);
        int comp1 = ComponentNumbering[first];
        int comp2 = ComponentNumbering[second];
        if (comp1 == comp2)
        {
        Graph::add_edge(components[comp1], correspondance[first], correspondance[second]);
        }
    }





//    std::vector<vertex_descriptor> v;
//    std::vector<std::pair<edge_descriptor, bool>> e;
//
//    size_t cloud_size = cloud.size();
//
//    for (int counter = 0; counter < cloud_size; ++counter)
//    {
//        v.push_back(Graph::add_vertex(component,cloud[counter].pt));
//   //     component[v[counter]].pt = cloud[counter].pt;
//    }
//
//    for (int iter_1 = 0; iter_1 < cloud_size; ++iter_1)
//    {
//        for (int iter_2 = iter_1 + 1; iter_2 < cloud_size; ++iter_2)
//        {
//            if (boost::edge( cloud[iter_1].index, cloud[iter_2].index, g ).second)
//            {
//                e.push_back(Graph::add_edge(component, iter_1, iter_2));
//               // boost::add_edge( iter_1, iter_2, component ) );
//                Point source = component[boost::source( e.back().first, component )].pt;
//                Point target = component[boost::target( e.back().first, component )].pt;
//                double length = norm( target - source );
//                boost::put( boost::edge_weight_t(), component, e.back().first, length );
//            }
//        }
//    }
}

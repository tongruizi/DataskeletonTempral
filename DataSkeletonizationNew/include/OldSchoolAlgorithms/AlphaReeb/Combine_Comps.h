#pragma once

#include "Graph.h"

void Combine_Comps ( std::vector<MyGraphType>const& conn_comp, MyGraphType& G )
{
    for (auto c_c : conn_comp)
    {
        std::vector<int> correspondance(boost::num_vertices(c_c));
        for (auto vi = boost::vertices( c_c ).first; vi != boost::vertices( c_c ).second; ++vi)
        {
            vertex_descriptor last = Graph::add_vertex(G,c_c[*vi].p);
            correspondance[*vi] = last;
        }
        for (auto eit = boost::edges(c_c).first; eit != boost::edges(c_c).second; eit++)
        {
            Graph::add_edge(G, correspondance[boost::source(*eit, c_c)],correspondance[boost::target(*eit, c_c)]);
        }
    }
}

//for (auto c_c : conn_comp)
//	{
//		std::vector<Graph::vertex_descriptor> v;
//		std::vector<std::pair<Graph::edge_descriptor, bool>> e;
//
//        size_t num_vertices = boost::num_vertices( c_c );
//
//		for (auto vi = boost::vertices( c_c ).first; vi != boost::vertices( c_c ).second; ++vi)
//		{
//			v.push_back( boost::add_vertex( g ) );
//			g[v.back()].pt = c_c[*vi].pt;
//		}
//
//        for (int counter_1 = 0; counter_1 < num_vertices; ++counter_1)
//        {
//            for (int counter_2 = counter_1 + 1; counter_2 < num_vertices; ++ counter_2)
//            {
//                if (boost::edge( counter_1, counter_2, c_c ).second)
//                {
//                    e.push_back( boost::add_edge( v[counter_1], v[counter_2], g ) );
//                    Point2d source = g[boost::source( e.back().first, g )].pt;
//                    Point2d target = g[boost::target( e.back().first, g )].pt;
//                    double length = norm( target - source );
//                    boost::put( boost::edge_weight_t(), g, e.back().first, length );
//                }
//
//            }
//        }
//	}

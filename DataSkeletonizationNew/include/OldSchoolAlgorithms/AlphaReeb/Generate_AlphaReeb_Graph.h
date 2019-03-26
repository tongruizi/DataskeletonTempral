#pragma once

#include "Graph.h"
#include "PointTools.h"

void Root ( std::vector<std::pair<bool, int>>const& tree, int& ptr )
{
    while (!tree[ptr].first)
    {
        ptr = tree[ptr].second;
    }
}

void Generate_AlphaReeb_Graph ( MyGraphType const& input_graph, double alpha, MyGraphType& alphaReeb_graph )
{
    size_t num_vertices = boost::num_vertices( input_graph );

    std::vector<std::pair<bool, int>> vertex_tree( 3 * num_vertices, std::pair<bool, int>( true, -1 ) );
    std::vector<std::pair<bool, int>> edge_tree( 2 * num_vertices, std::pair<bool, int>( true, -1 ) );

    for (auto ei = boost::edges( input_graph ).first; ei != boost::edges( input_graph ).second; ++ei)
    {
        int source, target;

        if (input_graph[boost::source( *ei, input_graph )].interval < input_graph[boost::source( *ei, input_graph )].interval)
        {
            source = boost::source( *ei, input_graph );
            target = boost::target( *ei, input_graph );
        }

        else
        {
            source = boost::target( *ei, input_graph );
            target = boost::source( *ei, input_graph );
        }

        int s = 2 * source + 1;
        int t = 2 * target;

        Root( edge_tree, s );
        Root( edge_tree, t );

        if (s != t)
        {
            edge_tree[t].first = false;
            edge_tree[t].second = s;
        }

        s = 3 * source + 1;
        t = 3 * target;

        Root( vertex_tree, s );
        Root( vertex_tree, t );

        if (s != t)
        {
            vertex_tree[t].first = false;
            vertex_tree[t].second = s;
        }

        s = 3 * source + 2;
        t = 3 * target + 1;

        Root( vertex_tree, s );
        Root( vertex_tree, t );

        if (s != t)
        {
            vertex_tree[t].first = false;
            vertex_tree[t].second = s;
        }
    }

    std::vector<vertex_descriptor> vd;

    for (int counter = 0; counter < 3 * num_vertices; ++counter)
    {
        if (vertex_tree[counter].first)
        {
            vertex_tree[counter].second = (int)vd.size();

            vd.push_back( Graph::add_vertex( alphaReeb_graph, Point(0,0,0) ) );

            if (counter % 3 == 1)
            {
                alphaReeb_graph[vd.back()].p = input_graph[(counter - 1) / 3].p;
            }

            else if (counter % 3 == 0)
            {
          //      Point qq = input_graph[counter / 3].p;
          //      Point qq2 = Point( 0, 0, (-1) * alpha / (double)2 );
           //      alphaReeb_graph[vd.back()].p = PointTools::sumPoint(qq,qq2);
           alphaReeb_graph[vd.back()].p = input_graph[counter / 3].p;
            }

            else
            {
           //     Point qq = input_graph[(counter - 2) / 3].p;
            //    Point qq2 = Point( 0, 0, alpha / (double) 2 );
            //    alphaReeb_graph[vd.back()].p = PointTools::sumPoint(qq,qq2);
             alphaReeb_graph[vd.back()].p = input_graph[(counter - 2) / 3].p;
            }
        }
    }

    for (int counter = 0; counter < 2 * num_vertices; ++counter)
    {
        if (edge_tree[counter].first)
        {
            int source, target;
            if (counter % 2 == 0)
            {
                source = counter * 1.5;
                target = counter * 1.5 + 1;
            }
            else
            {
                source = (counter - 1) * 1.5 + 1;
                target = (counter - 1) * 1.5 + 2;
            }
            Root( vertex_tree, source );
            Root( vertex_tree, target );
            source = vertex_tree[source].second;
            target = vertex_tree[target].second;
            Graph::add_edge(alphaReeb_graph , source, target);
        }
    }
}

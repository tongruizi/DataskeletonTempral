#pragma once
#include "Definitions.h"
#include<boost/graph/connected_components.hpp>
#include "ClusterElement.h"
#include "ClusteringMethods.h"

typedef std::pair<edge_descriptor,double> edgePair;

// is this really from biggest to smallest?

bool cmpedge(const edgePair &a, const edgePair &b)
{

    return a.second > b.second;
}

void SingleMSTAllocator(MyGraphType & mst, double param)
{
    std::vector<ClusterElement<edge_descriptor>> edges;
    boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, mst);

    for (auto eit = boost::edges(mst).first ; eit !=  boost::edges(mst).second; eit++)
    {
        edges.push_back(ClusterElement<edge_descriptor>(*eit, weightmap[*eit]));
    }

    ClusteringMethods<edge_descriptor> core;
    int j = core.DownwardClusterization(edges, param, "pure");

    if (j != -1)
    {
        for (int k = 0; k <= j; k++)
        {
        boost::remove_edge(edges[k].returnObject(),mst);
        }
    }


}


void AllocateComponents ( std::vector<MyGraphType>& msts, double param)
{
    for (int i = 0; i < msts.size(); i++)
    {
        SingleMSTAllocator(msts[i],param);
    }



    // First split MST into components:




    // Glue Components as it was done before



//    size_t num_intervals = subcloud.size();
//
//    vector<Graph> g(num_intervals);
//
//    Split_MST( subcloud, mst, g );
//
//    for (int counter_1 = 0; counter_1 < num_intervals; ++counter_1)
//    {
//        size_t num_vertices = boost::num_vertices( g[counter_1] );
//
//        vector<int> comp(num_vertices);
//        int num_comps = boost::connected_components( g[counter_1], &comp[0] ); // Assigns each vertex to its connected component.
//        vector<vector<Data_Pt>> conn_comp_cloud(num_comps);
//
//        for (int counter_2 = 0; counter_2 < num_vertices; ++counter_2)
//        {
//            conn_comp_cloud[comp[counter_2]].push_back( subcloud[counter_1][counter_2] );
//        }
//
//        for (int counter_2 = 0; counter_2 < num_comps; ++counter_2)
//        {
//            cluster.push_back( make_pair( conn_comp_cloud[counter_2], counter_1 ) );
//            cluster_vertex.push_back( Geometric_Centre_Of_Cloud( conn_comp_cloud[counter_2] ) );
//        }
//    }
}

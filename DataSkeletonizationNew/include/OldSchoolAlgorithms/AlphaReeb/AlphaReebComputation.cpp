#include "AlphaReebComputation.h"
#include "Dijkstra.h"
#include "Data_Pt.h"
#include "Split_Into_Conn_Comps.h"
#include "GraphUpgrade.h"
#include <map>
#include "Generate_AlphaReeb_Graph.h"
#include "Combine_Comps.h"
#include "GeneralConvertor.h"
#include "Graph.h"

AlphaReebComputation::AlphaReebComputation()
{
    //ctor
}

//int Generate_Subgraphs_Directly ( vector<Data_Pt>const& cloud, multimap<double, int>& filter_multimap, double alpha, vector<vector<Data_Pt>>& subcloud
//std::vector<int> & correspondance)

void writeDebug(std::string folder, std::vector<std::vector<Point>> & debugPointCloud)
{
    for (int i = 0; i < debugPointCloud.size(); i++)
    {
        MyGraphType G;
        for (auto it = debugPointCloud[i].begin(); it != debugPointCloud[i].end(); it++)
        {
            Graph::add_vertex(G,*it);
        }

        std::string path = folder + "debug" + std::to_string(i) + ".vtk";
        GeneralConvertor::GraphToVtk(path,G);
    }
}




void AlphaReebComputation::Compute( MyGraphType const& input_graph, AlphaReeb_Parameters const& parameters, MyGraphType& alphaReeb_graph,MyGraphType& InterMediate)
{
    int num_comps;
    std::vector<MyGraphType> conn_comp;
    Split_Into_Conn_Comps( input_graph, num_comps, conn_comp );
    std::vector<MyGraphType> alphaReeb_comp( num_comps );
    double min_comp_size = boost::num_vertices( input_graph ) * parameters.mcsf;
    for (int counter = 0; counter < num_comps; ++counter)
    {
        MyGraphType intermediate_graph;
        std::vector<int> intervalAllocation;
        int numberofvertices = boost::num_vertices(conn_comp[counter]);
        intervalAllocation.resize(numberofvertices);
        std::vector<std::vector<Point>> debugPointCloud;

        std::multimap<double, int> dijkstra_multimap; // Key = distance, value = index.
        Dijkstra( conn_comp[counter], dijkstra_multimap ); // Assigning filter values to each point.

        Generate_Subclouds_Correctly (conn_comp[counter], dijkstra_multimap, parameters.alpha, intermediate_graph, debugPointCloud);
        std::string debugfolder = "/home/yury/Dropbox/MicelleProject/Micelle/output/AlphaReeb/debug/";

     //   writeDebug(debugfolder, debugPointCloud);
        Generate_AlphaReeb_Graph( intermediate_graph, parameters.alpha, alphaReeb_comp[counter]);

        InterMediate = intermediate_graph;

    }
    Combine_Comps( alphaReeb_comp, alphaReeb_graph );

}




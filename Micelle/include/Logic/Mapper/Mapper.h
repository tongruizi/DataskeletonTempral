#pragma once

#include "Mapper_Parameters.h"
#include "Gauss_Density_Filter.h"
#include "Distance_Filter.h"
#include "Generate_Subclouds.h"
#include "Computation.h"
#include "ComponentAllocator.h"
#include "GraphUpgrade.h"
#include "Write.h"

double ComputeValue(std::vector<MyGraphType> trees, double scale)
{
    double countt = 0;
    double sum = 0;

    for (int i = 0; i < trees.size(); i++)
    {
        boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, trees[i]);
        for (auto eit = boost::edges(trees[i]).first; eit != boost::edges(trees[i]).second; eit++)
        {
            sum = sum + weightmap[*eit];
            countt++;
        }
    }
    return (sum/countt)*scale;
}

// void Computation::MSTSpecialCompute(std::map<Point, MyGraphType::vertex_descriptor> & vertex_map, MyGraphType & savedtree, std::vector<Point> & Vector)

void GenerateMST(std::vector<std::vector<Point>> & subcloud, std::vector<MyGraphType> & graphs, std::vector<std::map<Point, MyGraphType::vertex_descriptor>> & vertex_map )
{

    for (int i = 0; i < subcloud.size(); i++)
    {
        MyGraphType savedtree;
        Computation::MSTSpecialCompute(vertex_map[i], savedtree, subcloud[i]);
        graphs.push_back(savedtree);
    }
}

void RepairMapping( std::vector<std::map<Point, MyGraphType::vertex_descriptor>> & vertex_map, std::vector<std::unordered_map<int,int>>& vertexIntervals, std::vector<Point>const& cloud)
{
    for (int i = 0; i < cloud.size(); i++)
    {
        for (auto it = vertexIntervals[i].begin(); it != vertexIntervals[i].end(); it++)
        {
            (*it).second = vertex_map[it -> first][cloud[i]];

        }
    }
}




void debugMethod(std::vector<MyGraphType> & graphs, std::string folder)
{
    for (int i = 0; i < graphs.size(); i++)
    {
        Write::GraphToVtk(folder + "debug" + std::to_string(i) + ".vtk",graphs[i]);
    }
}

void debugPrinting(std::vector<std::unordered_map<int,int>> & VertexIntervals, std::vector<MyGraphType> & Graphs, std::vector<Point> const & cloud, std::string directory, std::string directory2)
{
    std::ofstream mystream;
    mystream.open(directory);
    std::ofstream mystream2;
    mystream2.open(directory2);
    for(int i = 0; i < VertexIntervals.size(); i++)
    {
        mystream << "Vertex: " << i << " Point: " << cloud[i] << std::endl;
        for (auto it = VertexIntervals[i].begin(); it != VertexIntervals[i].end(); it++)
        {
            mystream << it -> first << std::endl;
            mystream2 << "Compare: " << cloud[i] << " with: " << Graphs[it->first][it->second].p << std::endl;
        }
    }
    mystream.close();
}





void Mapper ( std::vector<Point>const& cloud, Mapper_Parameters const& parameters, MyGraphType& mapper_graph )
{
    // Assigning the filter value to each point.
    std::vector<std::map<Point, MyGraphType::vertex_descriptor>> vertex_map(parameters.num_intervals);
    std::vector<std::unordered_map<int,int>> vertexIntervals(cloud.size());

    std::multimap<double, int> filter_multimap; // Key = Gaussian density, value = index.
    filter_multimap.clear();

    if (parameters.filter_function == "Gauss_Density")
        Gauss_Density_Filter( cloud, parameters.sigma, filter_multimap );

    else if (parameters.filter_function == "Distance")
        Distance_Filter( cloud, filter_multimap );

    // Assigning points to subclouds.

    std::vector<std::vector<Point>> subcloud;
    subcloud.clear();
    subcloud.resize( parameters.num_intervals );

    Generate_Subclouds( cloud, filter_multimap, parameters.num_intervals, parameters.overlap_ratio, subcloud, vertexIntervals );

    // void Form_the_Graph(MyGraphType & G, MyGraphType & out, std::vector<MyGraphType> & Graphs, std::vector<std::unordered_map<int,int>> & vertexIntervals)

    // Creating the MST for each subcloud.

    std::vector<MyGraphType> mst;
    GenerateMST(subcloud,mst,vertex_map);
    double vv = ComputeValue(mst,parameters.mcsf);
    AllocateComponents ( mst, vv);
     RepairMapping(vertex_map,  vertexIntervals, cloud);

    // debug:
  //  std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/output/Mapper/debug/";
  // debugMethod(mst,folder);

  //  debugPrinting(vertexIntervals, mst, cloud, "/home/yury/Dropbox/MicelleProject/Micelle/output/Mapper/ddprint.txt",
           //       "/home/yury/Dropbox/MicelleProject/Micelle/output/Mapper/pointtest.txt");
    Form_the_Graph(cloud.size(), mapper_graph, mst, vertexIntervals);


}

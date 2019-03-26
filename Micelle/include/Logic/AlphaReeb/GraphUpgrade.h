#ifndef GraphUpgrade_H
#define GraphUpgrade_H

#include "PointTools.h"
#include "Graph.h"
#include <unordered_map>

// Alpha-Reeb.

void Intervals ( double alpha, double max, std::vector<std::pair<double, double>>& intervals );
void Form_the_Graph(int numVertices, MyGraphType & out, std::vector<MyGraphType> & Graphs, std::vector<std::unordered_map<int,int>> & vertexIntervals);
void Generate_Subclouds_Correctly (MyGraphType & G, std::multimap<double, int>& filter_multimap, double alpha, MyGraphType & out,
                                   std::vector<std::vector<Point>> & debugPointCloud);


#endif

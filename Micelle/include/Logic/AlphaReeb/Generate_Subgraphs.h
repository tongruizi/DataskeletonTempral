#pragma once

#include "Recover_Components.h"

void Generate_Subgraphs ( MyGraphTYpe const& g, std::vector<int> const& correspondance, vector<MyGraphType>& subgraph, int k )
{
    subgraph.resize(k);
	Recover_Components( g, correspondance, subgraph);

}

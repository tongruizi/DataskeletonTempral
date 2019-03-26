#pragma once

#include <boost/graph/connected_components.hpp>

#include "Recover_Components.h"

void Split_Into_Conn_Comps ( MyGraphType const& g, int& num_comps, std::vector<MyGraphType>& conn_comp )
{
	size_t num_vertices = boost::num_vertices( g );

	std::vector<int> comp( num_vertices );


	num_comps = boost::connected_components( g, &comp[0] ); // Assigns each vertex to its connected component.

		conn_comp.resize( num_comps );

 //   for (auto vi = boost::vertices( g ).first; vi != boost::vertices( g ).second; ++vi)
  //  {
   //     Data_Pt data_pt( g[*vi].pt );
    //    data_pt.index = g[*vi].index;
    //    conn_comp_cloud[comp[data_pt.index]].push_back( data_pt );
  //  }


    Recover_Components( g,  comp, conn_comp); // Assigns an edge between two vertices in each component if an edge exists in the input graph.

}

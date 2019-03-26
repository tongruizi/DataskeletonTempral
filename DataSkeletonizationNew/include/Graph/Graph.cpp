#include "Graph.h"
#include "Definitions.h"
#include "Computation.h"


Graph::Graph()
{
}
//  Point p;
//  double potentialvalue;
//  double maximumlength;
//  bool branchedcandidate;

vertex_descriptor add_vertexPrivate(MyGraphType &G, Point p)
{
vertex_descriptor v = add_vertex(G);
G[v].p = p;
G[v].index = v;
return v;
}

vertex_descriptor Graph::add_vertex(MyGraphType &G, Point p)
{
    return add_vertexPrivate(G,p);
}


edge_descriptor add_edgePrivate(MyGraphType &G,
                         vertex_descriptor v1,
                         vertex_descriptor v2, double d)
{
    edge_descriptor e = add_edge(v1, v2, G).first;
    boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, G);
    weightmap[e] = d;
    G[e].distance = d;
    return e;
}

edge_descriptor Graph::add_edge(MyGraphType &G,
                         vertex_descriptor v1,
                         vertex_descriptor v2)
{
double d = Computation::distance(G,v1,v2);
return add_edgePrivate(G,v1,v2,d);
}

edge_descriptor Graph::add_custom_edge(MyGraphType &G,
                         vertex_descriptor v1,
                         vertex_descriptor v2, double d)
{
return add_edgePrivate(G,v1,v2,d);
}


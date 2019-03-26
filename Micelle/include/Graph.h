#ifndef GRAPH_H
#define GRAPH_H
#include"Definitions.h"

struct VertexData;
struct EdgeData;

typedef boost::adjacency_list<boost::hash_setS, boost::vecS,
        boost::undirectedS,
        VertexData,
        boost::property<boost::edge_weight_t, double, EdgeData>
        > MyGraphType;

typedef typename boost::graph_traits<MyGraphType>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<MyGraphType>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<MyGraphType>::vertex_iterator vertex_iter;

struct EdgeData
{
    double distance;
};
struct VertexData
{
    Point p;
    vertex_descriptor index;
    vertex_descriptor correspondance;
    double potentialvalue;
    double maximumlength;
    bool branchedcandidate;
    bool boundary;
    int interval;
    vertex_descriptor prev;
    vertex_descriptor prev2;
    std::vector<std::pair<vertex_descriptor,MyGraphType*>> connectors;

};

class Graph
{
public:
    Graph();
    static vertex_descriptor add_vertex(MyGraphType &G, Point p);
    static edge_descriptor add_edge(MyGraphType &G,
                         vertex_descriptor v1,
                        vertex_descriptor v2);
    static edge_descriptor add_custom_edge(MyGraphType &G, vertex_descriptor v1,
    vertex_descriptor v2, double d);



protected:

private:
};

#endif // GRAPH_H

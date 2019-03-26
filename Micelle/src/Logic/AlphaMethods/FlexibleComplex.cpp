#include "FlexibleComplex.h"
#include "AABBMethods.h"
#include "Write.h"
#include "HomologyComputator.h"
#include "PointInfo.h"
#include "Computation.h"

FlexibleComplex::FlexibleComplex()
{

}

//void FlexibleComplex::ComputeCriticalGraph(MyGraphType & G, std::map<Point, PointInfo> & pointDistance)
//{
//
//    auto epair = boost::edges(savedTree);
//
//
//    unordered_set<int> criticalIndices;
//    std::list<std::pair<int,int>> Edges;
//    unordered_set<int, int> Correspondance;
//
//    auto vpair = boost::edges(savedTree);
//    for (auto it = vpair.first; it != vpair.second ; it++)
//    {
//    vertex_descriptor v1 = boost::source(*it, G);
//    vertex_descriptor v2 = boost::target(*it, G);
//    Point p1 = G[v1].p;
//    Point p2 = G[v2].p;
//    double abd = abs(pointDistance[p1] - pointDistance[p2]);
//
//    }
//
//}


void FlexibleComplex::MakePointDistanceMap(Alpha_shape_3 & as, std::map<Point, PointInfo> & pointDistance)
{
    std::list<Point> nonBoundary;
    std::list<Triangle> triangleList;
    int indexx = 0;
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        Point p = v_it -> point();
        pointDistance[p] = PointInfo{0,indexx};
        if (as.classify(v_it) == Alpha_shape_3::INTERIOR)
        {
            nonBoundary.push_back(p);
        }
        indexx++;
    }

    for (auto face = as.finite_facets_begin(); face != as.finite_facets_end(); face++)
    {
        if (as.classify(*face) == Alpha_shape_3::REGULAR)
        {
            auto tmp_tetra = face->first;
            Point apm[3];
            int j = 0;
            for (int i = 0; i < 4; i++)
            {
                if (i != face->second)
                {
                    apm[j] = tmp_tetra->vertex(i)->point();
                    j++;
                }
            }
            Triangle k(apm[0],apm[1],apm[2]);
            triangleList.push_back(k);
        }

    }
    AABBMethods aabb(triangleList);
    aabb.ComputeDistances(nonBoundary, pointDistance);

}

MyGraphType FlexibleComplex::OptimizedSpanningTree(Alpha_shape_3 & as, std::string filename)
{
    std::map<Point, PointInfo> pointDistance;
    FlexibleComplex::MakePointDistanceMap(as, pointDistance);
    double maxdistance = 0;
    // iterate overmap;
    MyGraphType G;
    for (auto viter =  as.finite_vertices_begin(); viter != as.finite_vertices_end(); viter++)
    {
        Point p =  viter -> point();
        vertex_descriptor pumpum = Graph::add_vertex(G, p);
        pointDistance[p].index = pumpum;
        double locald = pointDistance[p].distance;
        maxdistance = std::max(maxdistance, locald);
    }
    maxdistance = maxdistance + 1;
    for(auto edge = as.finite_edges_begin();
            edge != as.finite_edges_end();
            edge++)
    {
        bool k = as.classify(*edge) == Alpha_shape_3::EXTERIOR;
        if (k == false)
        {
            auto tmp_tetra = (*edge).get<0>();
            int p1, p2;
            p1 = (*edge).get<1>();
            p2 = (*edge).get<2>();
            PointInfo kk1 = pointDistance[tmp_tetra->vertex(p1)->point()];
            PointInfo kk2 = pointDistance[tmp_tetra->vertex(p2)->point()];
            int v0 = kk1.index;
            int v1 = kk2.index;
            double distancez = maxdistance - (kk1.distance + kk2.distance)/2;
            Graph::add_custom_edge(G,v0,v1,distancez);
        }

    }
    std::list<boost::graph_traits<MyGraphType>::edge_descriptor> mst_kruskal;
    boost::kruskal_minimum_spanning_tree(G, std::back_inserter(mst_kruskal));
    MyGraphType savedTree;
    Computation::treeGraph(G, mst_kruskal, savedTree);


    // Rebalance edges:
    auto epair = boost::edges(savedTree);
    boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, savedTree);
    for (auto it = epair.first; it != epair.second; it++)
    {
        double oldd = savedTree[(*it)].distance;
        oldd = maxdistance - oldd;
        auto v1 = boost::source((*it),savedTree);
        auto v2 = boost::target((*it),savedTree);
        double e1 = sqrt(CGAL::squared_distance(savedTree[v1].p,savedTree[v2].p));
        oldd = oldd * e1;
        weightmap[(*it)] = e1;
        savedTree[(*it)].distance = e1;
    }

    return savedTree;
    //Write::GraphToVtk(filename,savedtree);
}

void FlexibleComplex::AlphaTest(Alpha_shape_3 & as, std::string filename)
{
    int constant = 10;
    std::map<Point, PointInfo> pointDistance;

    Write::DistancesFromBoundary(filename,pointDistance);

    Simplex_tree st;
    // Adding all the simplices

    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        Point p = v_it -> point();
        PointInfo pi = pointDistance[p];
        auto apm = {pi.index};
        std::cout << "Vertex with distance: " << p << " : " << pi.distance << std::endl;
        st.insert_simplex(apm, constant-1*pi.distance);
    }
    for (auto edge = as.finite_edges_begin(); edge != as.finite_edges_end(); edge++)
    {
        bool k = as.classify(*edge) == Alpha_shape_3::EXTERIOR;
        if (k == false)
        {
            auto tmp_tetra = (*edge).get<0>();
            int p1, p2;
            p1 = (*edge).get<1>();
            p2 = (*edge).get<2>();
            PointInfo kk1 = pointDistance[tmp_tetra->vertex(p1)->point()];
            PointInfo kk2 = pointDistance[tmp_tetra->vertex(p2)->point()];
            int v0 = kk1.index;
            int v1 = kk2.index;
            double distancez = std::max(-1*kk1.distance, -1*kk2.distance);
            Chandler apm = {v0,v1};
            st.insert_simplex(apm, constant + distancez);
        }
    }
    for (auto face = as.finite_facets_begin(); face != as.finite_facets_end(); face++)
    {
        bool k = as.classify(*face) == Alpha_shape_3::EXTERIOR;
        if (k == false)
        {
            double d = std::numeric_limits<double>::min();
            auto tmp_tetra = face->first;
            int apm[3];
            int j = 0;
            for (int i = 0; i < 4; i++)
            {
                if (i != face->second)
                {
                    PointInfo tmpp = pointDistance[tmp_tetra->vertex(i)->point()];
                    apm[j] = tmpp.index;
                    d = std::max(d,-1*tmpp.distance);
                    j++;
                }
            }
            st.insert_simplex(apm, constant + d);
        }
    }
    for (auto cell = as.finite_cells_begin(); cell != as.finite_cells_end(); cell++)
    {
        if (as.classify(cell) != Alpha_shape_3::EXTERIOR)
        {
            double d = std::numeric_limits<double>::min();
            int apm[4];
            for (int i = 0; i < 4; i++)
            {
                PointInfo tmpp = pointDistance[cell->vertex(i)->point()];
                apm[i] = tmpp.index;
                d = std::max(d,-1*tmpp.distance);

            }
            st.insert_simplex(apm, constant + d);
//             size_t v0 = points.find(cell->vertex(0)->point())->second;

        }
    }


    HomologyComputator::Compute(st);

}

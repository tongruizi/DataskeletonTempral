#include "BranchSimplification.h"
#include "StraighteningMethods.h"
#include "PointTools.h"

BranchSimplification::BranchSimplification()
{
    //ctor
}

vertex_descriptor FindPointOperation(Point p, std::map<Point, vertex_descriptor> & PointToVertex, MyGraphType & G)
{
    auto ft = PointToVertex.find(p);
    if (ft == PointToVertex.end())
    {
        vertex_descriptor v = Graph::add_vertex(G,p);
        PointToVertex[p] = v;
        return v;
    }
    else
    {
        return (*ft).second;
    }
}

double ComputeDistance(std::list<Point> & s)
{
    Point prev;
    double totalDistance = 0;
    for (auto it = s.begin(); it != s.end(); it++)
    {
        if (it != s.begin())
        {
            totalDistance = totalDistance + sqrt(CGAL::squared_distance(*it,prev));
        }
        prev = *it;
    }
    return totalDistance;

}
//void BranchDetection::SimplifyBranchPoints(MyGraphType & G, std::list<std::pair<Point, std::list<Point>>> & outputInfo)
//{
//std::cout << "SimplifyBranchPoint started" << std::endl;
//    auto data = branchIndexation("",G);
//// Lets Go thought roots xD
//std::cout << "BranchDetection succeful" << std::endl;
//
//    for (auto it = (*(std::get<2>(data))).begin(); it != (*(std::get<2>(data))).end(); it++)
//    {
//        double t = G[*it].potentialvalue/(G[*it].potentialvalue + G[G[*it].prev2].potentialvalue);
//        Point s = PointTools::convexCombo(G[*it].p, G[G[*it].prev2].p,t);
//        std::list<Point> pointlist;
//        BranchDetection::PointCollector(*it,std::get<0>(data),pointlist);
//        outputInfo.push_back(std::make_pair(s,pointlist));
//    }
//// Return point:
//std::cout << "SimplifyBranchPoint ended" << std::endl;
//}

void BranchSimplification::Allocation(MyGraphType & G, std::list<std::pair<Point, std::list<Point>>> & outputInfo)
{
    std::vector<int> componentMap(num_vertices(G));
    int num = boost::connected_components(G, componentMap.data());
    std::list<Point> pointLists[num];
    auto vpair = boost::vertices(G);
    for(auto it = vpair.first; it != vpair.second; it++)
    {
        int apm = componentMap[*it];
        pointLists[apm].push_back(G[*it].p);
//   std::cout << "Adding " << G[*it].p << " To: "<< apm << std::endl;

    }
    for (int i = 0; i < num; i++)
    {
        Point q = PointTools::findBarycenter(pointLists[i]);
//    std::cout << "Component " << i << " Size: " << pointLists[i].size() << " Barycenter: " << q << std::endl;

        outputInfo.push_back(std::make_pair(q, pointLists[i]));
    }

}

void BranchSimplification::PathToGraph(MyGraphType & G, std::list<std::list<Point>> & in)
{
    std::map<Point, vertex_descriptor> PointToVertex;
    for(auto it = in.begin(); it != in.end(); it++)
    {
        vertex_descriptor down =  FindPointOperation((*it).front(), PointToVertex,G);
        vertex_descriptor up   =  FindPointOperation((*it).back(),  PointToVertex,G);
        double dd = ComputeDistance(*it);
        Graph::add_custom_edge(G,down,up,dd);
    }



}


void BranchSimplification::PathToGraphProper(MyGraphType & G, std::list<std::list<Point>> & in)
{
    std::map<Point, vertex_descriptor> PointToVertex;
    for(auto it = in.begin(); it != in.end(); it++)
    {
        vertex_descriptor core;
        for (auto git = (*it).begin(); git != (*it).end(); git++)
        {
            if (git != (*it).begin())
            {
                vertex_descriptor cur = FindPointOperation(*git, PointToVertex,G);
                Graph::add_edge(G, core, cur);
                core = cur;
            }
            else
            {
                core = FindPointOperation(*git, PointToVertex,G);
            }
        }
        vertex_descriptor down =  FindPointOperation((*it).front(), PointToVertex,G);
        vertex_descriptor up   =  FindPointOperation((*it).back(),  PointToVertex,G);
        double dd = ComputeDistance(*it);
        Graph::add_custom_edge(G,down,up,dd);
    }

}

void BranchSimplification::SimplifyIt(std::list<std::list<Point>> & in, std::list<std::list<Point>> & out, double epsilon)
{


    std::map<Point, vertex_descriptor> PointToVertex;
    MyGraphType G;
    for(auto it = in.begin(); it != in.end(); it++)
    {
        vertex_descriptor down =  FindPointOperation((*it).front(), PointToVertex,G);
        vertex_descriptor up   =  FindPointOperation((*it).back(),  PointToVertex,G);
        double dd = ComputeDistance(*it);

        if (dd < epsilon)
        {
            Graph::add_custom_edge(G,down,up,dd);
        }
        else
        {
        }
    }
    std::list<std::pair<Point, std::list<Point>>> outputInfo;
    BranchSimplification::Allocation(G, outputInfo);
    for (auto it = outputInfo.begin(); it != outputInfo.end(); it++)
    {
//        vertex_descriptor v = PointToVertex[(*it).first];
        for(auto bt = (*it).second.begin(); bt != (*it).second.end(); bt++)
        {
            // Something must be written in this loop;
            if (PointToVertex.find(*bt) == PointToVertex.end())
            {
            }
            vertex_descriptor u = PointToVertex[(*bt)];
            G[u].p = (*it).first;

        }
    }
    std::cout << "IN INFORMATION: " << std::endl;
    for(auto it = in.begin(); it != in.end(); it++)
    {
        std::cout << "Size: " << (*it).size() << std::endl;
        Point p1 = G[PointToVertex[(*it).front()]].p;
        Point p2 = G[PointToVertex[(*it).back()]].p;
        if (p1 != p2)
        {
            (*it).pop_front();
            (*it).push_front(p1);
            (*it).pop_back();
            (*it).push_back(p2);
            out.push_back(*it);

        }
    }
    std::cout << "OUT INFORMATION: " << std::endl;
    for (auto it = out.begin(); it != out.end(); it++)
    {

        std::cout << "Size: " << (*it).size() << std::endl;


    }


}

vertex_descriptor pickTheOtherOne(MyGraphType & G, vertex_descriptor k, edge_descriptor e)
{
    vertex_descriptor one = boost::source(e, G);
    vertex_descriptor two = boost::target(e, G);
    if (one == k)
    {
        return two;
    }
    else
    {
        return one;
    }

}


edge_descriptor catchTheOnlyEdge(MyGraphType & G, vertex_descriptor & iterating)
{
    for (auto eit = boost::out_edges(iterating, G).first ; eit !=  boost::out_edges(iterating, G).second; eit++)
    {
        return *eit;
    }
}

void AllocateByMST(std::vector<vertex_descriptor> & branchpoints, std::vector<bool> & branchIndicator, std::vector<std::list<Point>> & pointclouds, MyGraphType G)
{


    for (auto it = branchpoints.begin(); it != branchpoints.end(); it++)
    {
        std::deque<vertex_descriptor> descriptors;
//  std::list<edge_descriptor> edgesToBeRemoved;
        for (auto eit = boost::out_edges(*it, G).first; eit != boost::out_edges(*it, G).second; eit++)
        {
            vertex_descriptor other =  pickTheOtherOne(G, *it, *eit);
            if (branchIndicator[other] == 0)
            {
                descriptors.push_back(other);
            }
            // edgesToBeRemoved.push_back(*eit);
            boost::clear_vertex(*it, G);
        }
        while (!descriptors.empty())
        {
            vertex_descriptor iterating = descriptors.front();
            descriptors.pop_front();
            while(boost::out_degree(iterating, G) == 1)
            {
                pointclouds[*it].push_back(G[iterating].p);
                edge_descriptor catchedEdge = catchTheOnlyEdge(G, iterating);
                boost::remove_edge(catchedEdge,G);
                vertex_descriptor nexttVertex = pickTheOtherOne(G, iterating, catchedEdge);
                iterating = nexttVertex;
            }
            pointclouds[*it].push_back(G[iterating].p);

            for (auto vit = boost::adjacent_vertices(iterating, G).first ; vit != boost::adjacent_vertices(iterating, G).second; vit++)
            {
                descriptors.push_back(*vit);
            }
            boost::clear_vertex(iterating,G);

        }

//    std::deque<vertex_descriptor> container;
//    for (auto it = boost::vertices(G).first ; it != boost::vertices(G).second; it++)
//    {
//        int out_degree = boost::out_degree(*it,G);
//        if (out_degree == 1)
//        {
//            container.push_back(*it);
//        }
//    }
//    while (!container.empty())
//    {
//    vertex_descriptor iterating = container.front();
//    container.pop_front();
//    std::list<Point>
//    while(boost::out_degree(iterating, G) == 1)
//    {
//
//    }
//    }


    }
}




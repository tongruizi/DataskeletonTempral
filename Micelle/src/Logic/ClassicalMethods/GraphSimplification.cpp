#include "GraphSimplification.h"

GraphSimplification::GraphSimplification()
{
    //ctor
}

bool GraphSimplification::CorrectNumberOfEnds(int k, MyGraphType & G)
{
    int endcount = 0;

    for (auto it = boost::vertices(G).first; it != boost::vertices(G).second; it++)
    {
        int degree = boost::degree(*it,G);
        if (degree == 1)
        {
            endcount++;
        }
    }
    return (endcount == k);
}

bool GraphSimplification::RecognizeStraGraph(MyGraphType & G, int n)
{
//std::cout << "Parameter N:" << n << std::endl;
    auto vpair = boost::vertices(G);
    int degree1 = 0;
    int starn = 0;
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        int tmpD = boost::degree(*it, G);
        if (tmpD == 1)
        {
            degree1++;
        }
        else if (tmpD == n)
        {
            starn++;
        }
        else if ((tmpD != 2) && (tmpD != n) && (tmpD != 1))
        {
          //  std::cout << "Unsucceful, tmpD =" << tmpD << std::endl;
            return false;
        }
    }
 //   std::cout << "Stats: N-star: " << starn << " | 1-vertex: " << degree1 << std::endl;
    return (starn == 1) && (degree1 == n);
}

bool GraphSimplification::RecognizeDoubleStarGraph(MyGraphType & G, int n, int m)
{
std::set<int> setInt;
    int n2 = n + 1;
    int m2 = m + 1;
    setInt.insert(n2);
    setInt.insert(m2);
    int onev;
    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        int tmpD = boost::degree(*it, G);
        if (tmpD > 2)
        {
            auto bt = setInt.find(tmpD);
            if (bt == setInt.end())
            {
            std::cout << tmpD << std::endl;
                return false;
            }
            else
            {
      //      std::cout << "ERASED SOMETHING" << std::endl;
            setInt.erase(tmpD);
            }
        }
        else if (tmpD == 1)
        {
            onev++;

        }

    }
 //   std::cout << "SETINT :" << setInt.size() << " | " << "ONEV:" << onev << std::endl;
    return (setInt.size() == 0) && (onev == n + m);

}

vertex_descriptor pickTheOtherOne2(MyGraphType & G, vertex_descriptor k, edge_descriptor e)
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


edge_descriptor catchTheOnlyEdge2(MyGraphType & G, vertex_descriptor & iterating)
{
    for (auto eit = boost::out_edges(iterating, G).first ; eit !=  boost::out_edges(iterating, G).second; eit++)
    {
        return *eit;
    }
}

void GraphSimplification::SimplifyGraph(MyGraphType G, MyGraphType & S)
{
    std::vector<bool> branchPoints(boost::num_vertices(G));
    std::vector<int> correspondance(boost::num_vertices(G));
    std::deque<vertex_descriptor> container;
    for (auto it = boost::vertices(G).first ; it != boost::vertices(G).second; it++)
    {
        int out_degree = boost::degree(*it,G);
        if (out_degree == 1)
        {
            container.push_back(*it);
            correspondance[*it] = -1;
        }
        if (out_degree != 2)
        {
            branchPoints[*it] = true;
            Point pp = G[*it].p;
            vertex_descriptor AA = Graph::add_vertex(S, pp);
            correspondance[*it] = AA;
        }
        else
        {
            branchPoints[*it] = false;
            correspondance[*it] = -1;
        }

    }


    while (!container.empty())
    {
        vertex_descriptor beginn = container.front();
        vertex_descriptor beginNewVertex ;
        if (correspondance[beginn] == -1)
        {
            beginNewVertex = Graph::add_vertex(S, G[beginn].p);
            correspondance[beginn] = beginNewVertex;
        }
        else
        {
            beginNewVertex = correspondance[beginn];
        }

        container.pop_front();
        vertex_descriptor nextt;
        if (boost::out_degree(beginn,G) == 0)
        {
            continue;
        }
        vertex_descriptor prev = beginn;
        edge_descriptor grabber;
        for (auto eit = boost::out_edges(beginn, G).first ; eit !=  boost::out_edges(beginn, G).second; eit++)
        {
            grabber = *eit;
            nextt = pickTheOtherOne2(G, beginn,*eit);
        }

        boost::remove_edge(grabber, G);

        // We continue next after next
        while(branchPoints[nextt] == false)
        {
            edge_descriptor savedone;
            for (auto eit = boost::out_edges(nextt,G).first ; eit != boost::out_edges(nextt,G).second; eit++)
            {
                vertex_descriptor cur = pickTheOtherOne2(G, nextt, *eit);
                savedone = *eit;
                nextt = cur;
            }
            boost::remove_edge(savedone,G);
        }

        vertex_descriptor pairdd = correspondance[nextt];
        Graph::add_edge(S,pairdd, beginNewVertex );

        if (boost::out_degree(nextt, G) == 1)
        {
            container.push_back(nextt);
        }


    }
}


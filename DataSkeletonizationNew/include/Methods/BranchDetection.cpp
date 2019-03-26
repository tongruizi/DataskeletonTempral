#include "BranchDetection.h"
#include "PointTools.h"
#include "ClusterElement.h"
#include "ClusteringMethods.h"

BranchDetection::BranchDetection()
{

}

bool Compare(VertexGraphPair & a, VertexGraphPair & b)
{
    MyGraphType* G = a.second;
    MyGraphType* K = b.second;
    return ((*G)[a.first]).potentialvalue > ((*K)[b.first]).potentialvalue;
}
bool cmpinternal(const VertexGraphPair &a, const VertexGraphPair &b)
{
    MyGraphType* G = a.second;
    MyGraphType* K = b.second;
    return ((*G)[a.first]).potentialvalue > ((*K)[b.first]).potentialvalue;

}
bool my_cmp2(const VertexGraphPair &a, const VertexGraphPair &b)
{
    MyGraphType* G = a.second;
    MyGraphType* K = b.second;
    vertex_descriptor a1 = ((*G)[a.first]).connectors[1].first;
    vertex_descriptor b1 = ((*K)[b.first]).connectors[1].first;
    return ((*G)[a1]).potentialvalue > ((*K)[b1]).potentialvalue;
}
void printResults(std::string directory, MyGraphType & G, std::vector<VertexGraphPair> & dataArray)
{
    std::ofstream textstream;
    textstream.open(directory);
    for (auto finit = dataArray.begin(); finit != dataArray.end(); finit++)
    {
        textstream << G[(*finit).first].p << ": " ;
        for (auto itm = G[(*finit).first].connectors.begin();
                itm !=  G[(*finit).first].connectors.end(); itm++)
        {
            textstream << G[(*itm).first].potentialvalue << " ";

        }
        textstream << std::endl;
    }
}
edge_descriptor getOnlyEdge(vertex_descriptor v, MyGraphType* G)
{

    auto containerEdgeN = out_edges(v, *G);
    auto pointerEdgeN = containerEdgeN.first;
    edge_descriptor edgebetween = *pointerEdgeN;
    return edgebetween;


}

std::tuple<MyGraphType,std::vector<VertexGraphPair>,std::list<vertex_descriptor>*>
branchIndexation(std::string directory,MyGraphType G)
{
    std::list<vertex_descriptor>* root = new std::list<vertex_descriptor>();
    std::vector<VertexGraphPair> vertexContainer;
    auto vpair = vertices(G);
    for (auto iter = vpair.first; iter != vpair.second; iter++)
    {
        int sizedges = out_degree(*iter,G);
        G[*iter].boundary = false;
        G[*iter].branchedcandidate = false;
        if (sizedges == 1)
        {
            auto container = out_edges(*iter, G);
            auto pointer = container.first;
            edge_descriptor e = *pointer;
            double d = G[e].distance;
            G[*iter].potentialvalue = d;
            G[*iter].maximumlength = 0;
            G[*iter].boundary = true;
            VertexGraphPair element = std::make_pair(*iter, &G);
            vertexContainer.push_back(element);
        }
        else if (sizedges >= 3)
        {
            G[*iter].branchedcandidate = true;
        }
    }

    // std::cout << "Approaching heap" << std::endl;
    // Building heap:
    VertexHeap heap = VertexHeap(Compare,vertexContainer);
    // Running the program:
    while(!heap.empty())
    {
        VertexGraphPair cur = heap.top();
        heap.pop();
        vertex_descriptor it = cur.first;
        int tempdegree = out_degree(it,G);
        if (tempdegree == 0)
        {
            continue;
        }
        edge_descriptor edgebetween = getOnlyEdge(it, &G);
        auto container = adjacent_vertices(it, G);
        auto pointer = container.first;
        vertex_descriptor neighbor = *pointer;
        remove_edge(edgebetween, G);
        int degree = out_degree(neighbor,G);
        double pot = G[it].potentialvalue;
        if (degree != 0)
        {
            (G[neighbor].connectors).push_back(cur);
        }
        if (degree == 0)
        {
            //  G[neighbor].maximumlength2 = pot;
            G[neighbor].prev2 = it;
            (*root).push_back(neighbor);
            //	 VertexGraphPair
        }
        if (degree == 1)
        {
            G[neighbor].maximumlength = pot;
            edge_descriptor e = getOnlyEdge(neighbor, &G);
            double itsdistance = G[e].distance;
            double totaldistance = itsdistance + pot;
            G[neighbor].potentialvalue = totaldistance;
            G[neighbor].prev = it;
            VertexGraphPair paird = std::make_pair(neighbor, &G);
            heap.push(paird);
        }
        else if (degree >= 2)
        {

        }

    }
    //  std::cout << "After Heap" << std::endl;

    std::vector<VertexGraphPair> dataArray;
    // Getting the results into Array and sorting them in such way that highest low value first.
    auto vpairt = vertices(G);
    for (auto iter = vpairt.first; iter != vpairt.second; iter++)
    {
        if (G[*iter].branchedcandidate == true)
        {
            std::sort(G[*iter].connectors.begin(),G[*iter].connectors.end(), cmpinternal);
            dataArray.push_back(std::make_pair(*iter, &G));
        }
    }
    // sorting it:

    //  std::cout << "Before dataArray" << std::endl;

    //  std::cout << "DataArraySize: " << dataArray.size() << std::endl;
    std::sort(dataArray.begin(), dataArray.end(), my_cmp2);
    // Printing the results:

    // std::cout << "Before print results" << std::endl;

    printResults(directory, G, dataArray);

    //   std::cout << "After print results" << std::endl;

    // std::cout << " Roots found: " << calculationend.size() << std::endl;
    return std::make_tuple(G,dataArray,root);
}

bool acceptedCandidate(MyGraphType & G, vertex_descriptor c, double minvalue)
{
    if (G[c].connectors.size() < 2)
    {
        return false;
    }
    return G[G[c].connectors[1].first].potentialvalue > minvalue;
}

void searchForRoot(MyGraphType & G, MyGraphType & K,
                   vertex_descriptor g, vertex_descriptor k, std::deque<std::pair<vertex_descriptor,vertex_descriptor>> & Container,
                   double minvalue)
{
    vertex_descriptor kRoot = k;
    vertex_descriptor curE = g;
    while((G[curE].boundary == false)&&(!acceptedCandidate(G,curE,minvalue)))
    {
        curE = G[curE].connectors[0].first;
        Point p = G[curE].p;
        vertex_descriptor nw = Graph::add_vertex(K,p);
        Graph::add_edge(K,kRoot,nw);
        kRoot = nw;
    }
    if (G[curE].boundary == false)
    {
        for (int i = 0; i < G[curE].connectors.size(); i++)
        {
            vertex_descriptor tmp = G[curE].connectors[i].first;
            if (G[tmp].potentialvalue >= minvalue)
            {
                Point p = G[tmp].p;
                vertex_descriptor tmp13 = Graph::add_vertex(K,p);
                Graph::add_edge(K,tmp13,kRoot);
                Container.push_back(std::make_pair(tmp,tmp13));
            }
        }
    }



}

void GraphSimplification(MyGraphType & G, MyGraphType & K, double minvalue, vertex_descriptor root)
{
    std::deque<std::pair<vertex_descriptor,vertex_descriptor>> Container;
    vertex_descriptor bor = root;
    Point abba = G[root].p;
    vertex_descriptor nr = Graph::add_vertex(K,abba);
    // add to container the base things
    vertex_descriptor or2 = G[root].prev2;
    vertex_descriptor nr2 = Graph::add_vertex(K,G[or2].p);
    Graph::add_edge(K,nr,nr2);
    Container.push_back(std::make_pair(bor,nr));
    Container.push_back(std::make_pair(or2,nr2));
    while(!Container.empty())
    {
        auto vpair = Container.front();
        Container.pop_front();
        searchForRoot(G,K,vpair.first, vpair.second, Container, minvalue);
    }

}

double minvalueFinder(MyGraphType & G, double scale, std::string settings)
{
    std::vector<ClusterElement<vertex_descriptor>> examElements;
    for (auto vit = boost::vertices(G).first; vit != boost::vertices(G).second; vit++)
    {
        if (G[*vit].connectors.size() > 1)
        {
            examElements.push_back(ClusterElement<vertex_descriptor>(*vit,G[G[*vit].connectors[1].first].potentialvalue));
        }
    }
    ClusteringMethods<vertex_descriptor> engine;
    int indexx = engine.DownwardClusterization(examElements,scale,settings) + 1;
    engine.debugSet(examElements, "");
    if (indexx == -1)
    {
        return 2*examElements[0].returnValue();
    }
    else
    {
        return examElements[indexx].returnValue();

    }
}
//void BranchDetection::Convert(MyGraphType & G, )

void BranchDetection::SimplifyIt(MyGraphType & G, MyGraphType & K, double scale, std::string debug, std::string settings)
{
    auto data = branchIndexation(debug, G);
    double minvalue = 0;
    if (settings == "pure")
    {
        minvalue =  scale;
    }
    else
    {
        minvalue = minvalueFinder(std::get<0>(data), scale, settings);
      //  std::cout << "minvalue: " << minvalue << std::endl;

    }
    for (auto it = (*(std::get<2>(data))).begin(); it != (*(std::get<2>(data))).end(); it++)
    {
        GraphSimplification(std::get<0>(data),K,minvalue,*it);
    }

}

void BranchDetection::PointCollector(vertex_descriptor k, MyGraphType & G, std::list<Point> & pointlist)
{

    std::deque<std::pair<vertex_descriptor,vertex_descriptor>> vertexDeque;
    vertexDeque.push_back(std::make_pair(k,k));
    vertexDeque.push_back(std::make_pair(G[k].prev2,G[k].prev2));
    while(!vertexDeque.empty())
    {
        vertex_descriptor it = vertexDeque.front().first;
        vertex_descriptor bt = vertexDeque.front().second;
        vertexDeque.pop_front();
        if (!G[it].boundary)
        {
            for (int i = 0; i < G[it].connectors.size(); i++)
            {
                vertex_descriptor atm = G[it].connectors[i].first;
                if (atm != bt)
                {
                    vertexDeque.push_back(std::make_pair(atm,it));
                }
            }
        }
        else
        {
            pointlist.push_back(G[it].p);
        }
    }




}

void BranchDetection::SimplifyBranchPoints(MyGraphType & G, std::list<std::pair<Point, std::list<Point>>> & outputInfo)
{
    auto data = branchIndexation("",G);
// Lets Go thought roots xD

    for (auto it = (*(std::get<2>(data))).begin(); it != (*(std::get<2>(data))).end(); it++)
    {
        double t = G[*it].potentialvalue/(G[*it].potentialvalue + G[G[*it].prev2].potentialvalue);
        Point s = PointTools::convexCombo(G[*it].p, G[G[*it].prev2].p,t);
        std::list<Point> pointlist;
        BranchDetection::PointCollector(*it,std::get<0>(data),pointlist);
        outputInfo.push_back(std::make_pair(s,pointlist));
    }
// Return point:
}


#include "MetricGraphRec/MetricRec.hpp"
#include <mlpack/methods/range_search/range_search.hpp>
#include "GeneralConvertor.h"
#include <unordered_map>


using namespace mlpack::range;


void MetricRec::Run(std::list<Point> & cloudlist, MyGraphType & out)
{
    this->Labeling(cloudlist,out);
//1!

}

void MetricRec::VectorToMat(arma::mat & original, arma::mat & newOne, std::vector<size_t> & indicesSelected)
{


    newOne.resize(original.n_rows, indicesSelected.size());
    for(int i = 0; i < indicesSelected.size(); i++)
    {
        newOne.col(i) = original.col(indicesSelected[i]);
    }
}

void MetricRec::MedialPointComputator(std::vector<int> & indices, arma::mat & matrix, arma::vec & fProduct)
{
    fProduct.resize(matrix.n_rows);
    double multiplier = 1 / (double) indices.size();
    for (int i = 0; i < indices.size(); i++)
    {
        fProduct = fProduct + matrix.col(indices[i]);
    }
    fProduct = multiplier * fProduct;

}

void MetricRec::Labeling(std::list<Point> & cloudlist, MyGraphType & out)
{
//get y←BY(y,5r/3)\BY(y,r)
    //! Edgepoint and branchpoints have to be declared here:
    std::unordered_set<int> edgePoint;
    std::unordered_set<int> branchPoint1;
    std::unordered_set<int> branchPoint;


    std::vector<Point> vectorList(cloudlist.begin(),cloudlist.end());
    //! Other stuff:
    arma::mat cordata;
    GeneralConvertor::ListToMatTransposed(cloudlist,cordata);
    mlpack::range::RangeSearch<> a(cordata);

    std::vector<std::vector<size_t> > resultingNeighbors;
    std::vector<std::vector<double> > resultingDistances;

    double lowbound = this->r;
    double highbound = 5*this->r/3;
    mlpack::math::Range t(lowbound, highbound);
    std::cout << "Lowbound: " << lowbound << std::endl;
    std::cout << "Highbound: " << highbound << std::endl;
    a.Search(t, resultingNeighbors, resultingDistances);

    for (int i = 0; i < resultingNeighbors.size(); i++)
    {
        MyGraphType G;
        for (int k = 0; k<resultingNeighbors[i].size(); k++)
        {
            Graph::add_vertex(G, vectorList[resultingNeighbors[i][k]]);

        }
        //connected components of Rips-Vietoris graphR4r/3(Sy)
        arma::mat cordata1;
        // GeneralConvertor::ListToMatTransposed(resultingNeighbors[i],cordata1);
        this->VectorToMat(cordata, cordata1, resultingNeighbors[i]);
        mlpack::range::RangeSearch<> b(cordata1);

        std::vector<std::vector<size_t> > resultingNeighbors1;
        std::vector<std::vector<double> > resultingDistances1;
        double multiplier = 4/3;
        mlpack::math::Range p(0,multiplier*this->r);
        b.Search(p, resultingNeighbors1, resultingDistances1);
        //degr(y)←Number of connected components
        for (int k = 0; k<resultingNeighbors1.size(); k++)
        {
            for (int q = 0; q<resultingNeighbors1[k].size(); q++)
            {
                Graph::add_edge(G, k,resultingNeighbors1[k][q]);
            }
        }
        std::vector<int> componentMap(num_vertices(G));
        int num = boost::connected_components(G, componentMap.data());
        std::cout << "Point " << i << "gets num:" << num << std::endl;
        //first label points by dgr(y)
        if(num != 2)
        {
            branchPoint1.insert(i);
        }

    }
    //Label all points within distance 2tfrom a preliminary branch point as branch points

    std::unordered_map<int,double> debugmap;
    std::vector<std::vector<size_t> > resultingNeighbors2;
    std::vector<std::vector<double> > resultingDistances2;
    mlpack::math::Range q(0,2*this->t);
    a.Search(q, resultingNeighbors2, resultingDistances2);

    for (int i = 0; i < resultingNeighbors2.size(); i++)
    {
        if (branchPoint1.find(i) != branchPoint1.end())
        {
            for (int q = 0; q<resultingNeighbors2[i].size(); q++)
            {

                branchPoint.insert(resultingNeighbors2[i][q]);

                if (debugmap.find(resultingNeighbors2[i][q]) == debugmap.end())
                {
                    debugmap[resultingNeighbors2[i][q]] = resultingDistances2[i][q];
                }
                else
                {
                    double dd = debugmap[i];
                    if (resultingDistances2[i][q] < dd)
                    {
                        debugmap[resultingNeighbors2[i][q]] = resultingDistances2[i][q];

                    }
                }
            }
        }
    }

    for (int i = 0; i < cloudlist.size(); i++)
    {
        if (debugmap.find(i) == debugmap.end())
        {
            std::cout << "information for " << i << " not found " << std::endl;

        }
        std::cout << "Point: " << i << " : " << debugmap[i] << std::endl;
    }


    std::cout << "Cloudlist size: " << cloudlist.size() << std::endl;
    std::cout << "Branchpoint size: " << branchPoint.size() << std::endl;

    for (int i = 0;  i < cloudlist.size(); i++)
    {
        if (branchPoint.find(i) == branchPoint.end())
        {
            edgePoint.insert(i);
        }
    }

    std::cout << "Labeling finnishsed succefully " << std::endl;
    //! Continue in next method:
    std::cout << "Number of edgePoints:" << edgePoint.size() << std::endl;
    std::cout << "Number of branchpoints: " << branchPoint.size() << std::endl;



    //! Debug coloring:
    std::vector<double> coloring(cloudlist.size());
    for (int i = 0; i < cloudlist.size(); i++)
    {
        if (branchPoint1.find(i) == branchPoint1.end())
        {
            coloring[i] = 0;
        }
        else
        {
            coloring[i] = 1;
        }

    }
    arma::mat wCopy = cordata.t();
    std::string outt = "/home/liudi/Dropbox/LocalTests/OutputMetricGraphRec/debugcolor.csv";
    GeneralConvertor::MatInfoToFile(outt,wCopy,coloring);

    //! Second type of coloring:

    std::vector<double> coloring2(cloudlist.size());
    for (int i = 0; i < cloudlist.size(); i++)
    {
        if (branchPoint.find(i) == branchPoint.end())
        {
            coloring2[i] = 0;
        }
        else
        {
            coloring2[i] = 1;
        }

    }
    std::string outt2 = "/home/liudi/Dropbox/LocalTests/OutputMetricGraphRec/debugcolorActual.csv";
    GeneralConvertor::MatInfoToFile(outt2,wCopy,coloring2);


    this->ReconstructGraph(resultingNeighbors2,cordata,branchPoint,edgePoint,out);
}

void MetricRec::ProcedureForComputingComponents(int cloudnumber, std::unordered_set<int> & branchPoints,
        std::vector<std::vector<int>> & pointsToComponents, std::vector<std::vector<size_t>> &resultingNeighbors)
{

//! Compute the connected components of the Vietoris rips Graph:

    AbstractGraphType bGraph;
    std::vector<int> vMap(branchPoints.size());
    std::vector<int> inversevMap(cloudnumber);
    for (auto it = branchPoints.begin(); it != branchPoints.end(); it++)
    {
        int v = boost::add_vertex(bGraph);
        vMap[v] = *it;
        inversevMap[*it] = v;
    }

    for (auto it = branchPoints.begin(); it != branchPoints.end(); it++)
    {

        for (int j = 0; j < resultingNeighbors[*it].size(); j++)
        {
            if (branchPoints.find(resultingNeighbors[*it][j]) != branchPoints.end())
            {
                boost::add_edge(inversevMap[resultingNeighbors[*it][j]], inversevMap[*it], bGraph);
            }

        }
    }
    std::cout << "Graph has vertices: " << boost::num_vertices(bGraph) << std::endl;
    std::cout << "Graph has edges: " << boost::num_edges(bGraph) << std::endl;
    std::vector<int> componentMap(boost::num_vertices(bGraph));
    int num = boost::connected_components(bGraph, componentMap.data());

    std::cout << "Num components: " << num << std::endl;

    //! For each component compute the medial point

    pointsToComponents.resize(num);

    for (int i = 0; i < componentMap.size(); i++)
    {
        pointsToComponents[componentMap[i]].push_back(vMap[i]);
    }



}


void MetricRec::ReconstructGraph(std::vector<std::vector<size_t>> &resultingNeighbors, arma::mat & data, std::unordered_set<int> & branchPoints, std::unordered_set<int> & edgePoints, MyGraphType & output)
{
    std::cout << "Graph reconstruction started " << std::endl;
    //! Compute the connected components of the Vietoris rips Graph:
    int cloudnumber = data.n_cols;
    std::vector<std::vector<int>> pointsToComponentsVertices;
    std::cout << "Prodcure started for branchpoints" << std::endl;
    this->ProcedureForComputingComponents(cloudnumber,branchPoints,pointsToComponentsVertices,resultingNeighbors);
    // MedialPointComputator(std::vector<int> & indices, arma::mat & matrix, arma::vec & fProduct)
    for (int i = 0; i < pointsToComponentsVertices.size(); i++)
    {
        arma::vec p;
        MedialPointComputator(pointsToComponentsVertices[i],data,p);
        Graph::add_vertex(output, Point(p(0),p(1),p(2)));
    }

    //! Compute connected components of the Edges:
    std::vector<std::vector<int>> pointsToComponentsEdges;
    std::cout << "Procedure started for edge points " << std::endl;
    this->ProcedureForComputingComponents(cloudnumber,edgePoints,pointsToComponentsEdges,resultingNeighbors);

    std::cout << "Number of Edge clusters: " << pointsToComponentsEdges.size()<< std::endl;
    std::cout << "Number of Vertex clusters: " << pointsToComponentsVertices.size() << std::endl;

    std::vector<int> ComponentNumberEdge(cloudnumber);
    for (int i = 0; i < pointsToComponentsEdges.size(); i++)
    {
        for (int j = 0; j < pointsToComponentsEdges[i].size(); j++)
        {
            ComponentNumberEdge[pointsToComponentsEdges[i][j]] = i;

        }

    }


    std::vector<int> ComponentNumberVertex(cloudnumber);
    for (int i = 0; i < pointsToComponentsVertices.size(); i++)
    {
        for (int j = 0; j < pointsToComponentsVertices[i].size(); j++)
        {
            ComponentNumberVertex[pointsToComponentsVertices[i][j]] = i;

        }

    }

    std::cout <<  ComponentNumberVertex[0]<<"debug"<< std::endl;
    //! What do we do now?
    std::vector<std::unordered_set<int>> Collector(pointsToComponentsEdges.size());

    for (auto it = edgePoints.begin(); it != edgePoints.end(); it++)
    {
        for (int j = 0; j < resultingNeighbors[*it].size(); j++)
        {

            if (branchPoints.find(resultingNeighbors[*it][j]) != branchPoints.end())
            {
                int tmpp = resultingNeighbors[*it][j];
                Collector[ComponentNumberEdge[*it]].insert(ComponentNumberVertex[tmpp]);

            }

        }
    }


    for (int i = 0; i < Collector.size(); i++)
    {
        for (auto it = Collector[i].begin(); it != Collector[i].end(); it++)
        {
            for (auto bt = Collector[i].begin(); bt != Collector[i].end(); bt++)
            {
                if (*bt != *it)
                {
                    Graph::add_edge(output, *bt, *it);

                }
            }
        }

    }

}
















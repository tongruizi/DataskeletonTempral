
#include "GraphUpgrade.h"
#include "PointTools.h"
#include "Graph.h"
#include <unordered_map>

// Alpha-Reeb.

void Intervals ( double alpha, double max, std::vector<std::pair<double, double>>& intervals )
{
    double interval_start = 0, interval_end;

    for (int counter = 0; interval_start < max; ++counter)
    {
        interval_start = 0.5 * counter * alpha;
        interval_end = (0.5 * counter + 1) * alpha;
        intervals.push_back( std::pair<double, double>( interval_start, interval_end ) );
    }
}



// Replace G with std::vector

// Do something with vertexIntervals

void Form_the_Graph(int numVertices, MyGraphType & out, std::vector<MyGraphType> & Graphs, std::vector<std::unordered_map<int,int>> & vertexIntervals)
{
    int num_intervals = Graphs.size();

    std::vector<std::vector<std::vector<Point>>> componentPoints(num_intervals);
    std::vector<std::vector<int>> correspondance(num_intervals);
    std::vector<std::vector<int>> comp(num_intervals);
    std::vector<int> sizeComp(num_intervals);
    for(int i = 0; i < num_intervals; i++)
    {
        comp[i] = std::vector<int>(boost::num_vertices(Graphs[i]));
        sizeComp[i] = boost::connected_components(Graphs[i], &(comp[i][0]));
        componentPoints[i] = std::vector<std::vector<Point>>(sizeComp[i]);
        correspondance[i].resize(sizeComp[i]);
    }

    for(int i = 0; i < num_intervals; i++)
    {
        auto vpair = boost::vertices(Graphs[i]);
        for (auto vit = vpair.first; vit != vpair.second; vit++)
        {
        int compNumber = comp[i][*vit];
        componentPoints[i][compNumber].push_back(Graphs[i][*vit].p);
        // Graph will be added
        }
        for (int j = 0; j < sizeComp[i]; j++)
        {
        Point barycenter = PointTools::findBarycenter(componentPoints[i][j]);
        vertex_descriptor ne = Graph::add_vertex(out,barycenter);
        out[ne].interval = i;
        correspondance[i][j] = ne;
        }
    }


    // Handle the edges:

    for (int it = 0; it < numVertices; it++)
    {
        // Double loop is allright here as the container contains at most 3 objects (usually at most 2 objects)

        for (auto gitbit = vertexIntervals[it].begin(); gitbit != vertexIntervals[it].end(); gitbit++)
        {
            for (auto bitgit = vertexIntervals[it].begin(); bitgit != vertexIntervals[it].end(); bitgit++)
            {
                int indx1 = (*gitbit).first;
                int indx2 = (*bitgit).first;

                if (std::abs(indx1-indx2) == 1)
                {
                int descriptor1 = (*gitbit).second;
                int descriptor2 = (*bitgit).second;
                int comp1 = comp[indx1][descriptor1];
                int comp2 = comp[indx2][descriptor2];

                Graph::add_edge(out, correspondance[indx1][comp1], correspondance[indx2][comp2]);
                }
                else
                {

                }
            }
        }

    }


}


void Generate_Subclouds_Correctly (MyGraphType & G, std::multimap<double, int>& filter_multimap, double alpha, MyGraphType & out,
                                   std::vector<std::vector<Point>> & debugPointCloud)
{
    std::vector<std::pair<double, double>> intervals;

    auto it = filter_multimap.rbegin();
    double max_dist_from_root = it->first;
    Intervals( alpha, max_dist_from_root, intervals ); // Finding intervals.
    size_t num_intervals = intervals.size();
    debugPointCloud.resize(num_intervals);
    // We will have num_intervalls ammount of intervals.

    std::multimap<double, int>::iterator it_start, it_end;
    std::vector<std::pair<std::multimap<double, int>::iterator, std::multimap<double, int>::iterator>> pointers;

    for (int counter = 0; counter < num_intervals; ++counter)
    {
        it_start = filter_multimap.lower_bound( intervals[counter].first );
        it_end = filter_multimap.upper_bound( intervals[counter].second );
        pointers.push_back( std::make_pair( it_start, it_end ) );
    }

    // Assigning points to subclouds.

    std::vector<std::unordered_map<int,int>> vertexIntervals;
    vertexIntervals.resize(boost::num_vertices(G));
    std::vector<MyGraphType> Graphs;
    std::vector<std::vector<std::vector<Point>>> componentPoints(num_intervals);
    std::vector<std::vector<int>> correspondance(num_intervals);
    Graphs.resize(num_intervals);

    for (int counter = 0; counter < num_intervals; ++counter)
    {
        for (auto it = pointers[counter].first; it != pointers[counter].second; ++it)
        {
            Point QaQ = G[(*it).second].p;
            vertex_descriptor newGuy = Graph::add_vertex(Graphs[counter],QaQ);
        //    Graphs[counter][(*it).second].correspondance = *it;
            vertexIntervals[it->second].insert({counter,newGuy}); // LOOK at this
            debugPointCloud[counter].push_back(QaQ);
            //	subcloud[counter].push_back( cloud[it->second] );
        }
    }


    auto epair = boost::edges(G);
    for (auto eit = epair.first; eit != epair.second; eit++)
    {
        vertex_descriptor q1 = boost::source(*eit,G);
        vertex_descriptor q2 = boost::target(*eit,G);
        for (auto bit = vertexIntervals[q1].begin() ; bit != vertexIntervals[q1].end() ; bit++)
        {
            auto didFind = vertexIntervals[q2].find((*bit).first);
            if (didFind != vertexIntervals[q2].end())
            {
                Graph::add_edge(Graphs[(*bit).first], (*bit).second, didFind -> second);

            }
        }
    }

    std::vector<std::vector<int>> comp(num_intervals);
    std::vector<int> sizeComp(num_intervals);
    for(int i = 0; i < num_intervals; i++)
    {
        comp[i] = std::vector<int>(boost::num_vertices(Graphs[i]));
        sizeComp[i] = boost::connected_components(Graphs[i], &(comp[i][0]));
        componentPoints[i] = std::vector<std::vector<Point>>(sizeComp[i]);
        correspondance[i].resize(sizeComp[i]);
    }


    // We know now the connected components, we proceed into the next step:

    // Point Geometric_Centre_Of_Cloud ( vector<Data_Pt>const& cloud )


    // Handles the vertices:


    for(int i = 0; i < num_intervals; i++)
    {
        auto vpair = boost::vertices(Graphs[i]);
        for (auto vit = vpair.first; vit != vpair.second; vit++)
        {
        int compNumber = comp[i][*vit];
        componentPoints[i][compNumber].push_back(Graphs[i][*vit].p);
        // Graph will be added
        }
        for (int j = 0; j < sizeComp[i]; j++)
        {

        Point barycenter = PointTools::findBarycenter(componentPoints[i][j]);
        vertex_descriptor ne = Graph::add_vertex(out,barycenter);
        out[ne].interval = i;
        correspondance[i][j] = ne;
        }
    }


    // Handle the edges:


    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        // Double loop is allright here as the container contains at most 3 objects (usually at most 2 objects)

        for (auto gitbit = vertexIntervals[*it].begin(); gitbit != vertexIntervals[*it].end(); gitbit++)
        {
            for (auto bitgit = vertexIntervals[*it].begin(); bitgit != vertexIntervals[*it].end(); bitgit++)
            {
                int indx1 = (*gitbit).first;
                int indx2 = (*bitgit).first;
                int descriptor1 = (*gitbit).second;
                int descriptor2 = (*bitgit).second;
                int comp1 = comp[indx1][descriptor1];
                int comp2 = comp[indx2][descriptor2];

                if (std::abs(indx1-indx2) == 1)
                {
                Graph::add_edge(out, correspondance[indx1][comp1], correspondance[indx2][comp2]);
                }
                else
                {

                }
            }
        }

    }





}



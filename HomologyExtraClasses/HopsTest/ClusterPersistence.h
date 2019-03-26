#ifndef CLUSTERPERSISTENCE_H
#define CLUSTERPERSISTENCE_H

#include "CleverSkeletonization.h"
#include "DensityComputator.h"
#include "GeneralConvertor.h"
#include "AmstStat.h"

#include "GudhiSettings.h"

template<typename DensityType>
class ClusterPersistence
{

public:

    // typedef  mlpack::tree::KDTree<> ourTree;
    typedef  mlpack::tree::KDTree<mlpack::metric::EuclideanDistance, AmstStat, arma::mat> ourTree;
    typedef  DensityComputator<mlpack::metric::EuclideanDistance, arma::mat, DensityType, mlpack::tree::KDTree> ourDensityComputator;



    ClusterPersistence()
    {


    }

    virtual ~ClusterPersistence()
    {


    }



    double nv(double k, double q, double p, double x)
    {
        double a = (k-1)/(q-p);
        double b = 1-a*p;
        return a*x+b;

    }

    void Rescale(std::vector<double> & f, double t)
    {
        double maxx = 0;
        double minn = DBL_MAX;
        for (int i = 0; i < f.size(); i++)
        {
            maxx = std::max(f[i],maxx);
            minn = std::min(f[i], minn);
        }
        for (int i = 0; i < f.size(); i++)
        {
            f[i] = nv(t, minn,maxx, f[i]);
        }

    }
    // int & num, std::vector<int> & componentMap

    void ComputeClustering(std::vector<Point> & points, std::vector<size_t> & clusterCorrespondance, std::vector<size_t> & leaders, double fEpsilon, double clusterEpsilon, int & num, std::vector<size_t> & componentMap, std::string colorout)
    {
        //! Initilization
        std::vector<size_t> oldfromnew(points.size());
        arma::mat MatPoints;
        GeneralConvertor::VectorToMatTransposed(points,MatPoints);

        //AmstStat
        ourTree* ktree = new ourTree(MatPoints, oldfromnew);
        std::vector<double> f(MatPoints.n_cols);
        std::vector<int> vn(MatPoints.n_cols);
        //! Density computation
        ourDensityComputator calcf(ktree);

        //     void ComputeDensity(std::vector<double> & values, std::vector<int> & visitNumber, double epsilon)

        calcf.ComputeDensity(f, vn, fEpsilon);
        double t = 2;
        //! Rescaling
        Rescale(f,t);
        //! Compute new MST
        CleverSkeletonization<mlpack::metric::EuclideanDistance, arma::mat, mlpack::tree::KDTree> finalCalc(ktree,oldfromnew, false);

        //! Colorut:

        //! f,clusterEpsilon
        finalCalc.ComputeClusterization(clusterCorrespondance,f, clusterEpsilon,num,componentMap);

        //! Printout
        std::vector<double> doubleVec(clusterCorrespondance.begin(), clusterCorrespondance.end());
        arma::mat TransposedPoints = MatPoints.t();
        GeneralConvertor::MatInfoToFile(colorout,TransposedPoints,doubleVec);

        //! Deletion:
        ktree = NULL;
        delete ktree;

    }

    void DefineMapToMajorElement(int num, std::vector<size_t> & componentMap, std::vector<size_t> & clusterCorrespondance, std::vector<size_t> & resultingMap)
    {
        resultingMap.resize(num);
        for (int i = 0; i < clusterCorrespondance.size(); i++)
        {
            resultingMap[componentMap[i]] = clusterCorrespondance[i];
        }
        for (int i = 0; i < num ; i++)
        {
            std::cout << "Resulting Map( " << i << " ) = " << resultingMap[i] << std::endl;


        }

    }

    void ComputeDeluanayTriangulation( MyGraphType & G, std::vector<size_t> & clusterCorrespondance, std::vector<size_t> & resultingMap, std::vector<size_t> & componentMap,std::vector<Point> & Points, int & num)
    {
        // Triangulation T(Vector.begin(),Vector.end());
        //  MyGraphType Tmp;
        Triangulation T(Points.begin(),Points.end());

        //! Construct new Graph:
        for (int i = 0; i < num; i++)
        {
            Graph::add_vertex(G,Points[resultingMap[i]]);
        }

        std::map<Point, MyGraphType::vertex_descriptor> vertex_map;
        Triangulation::Finite_vertices_iterator viter;
        Triangulation::size_type n = T.number_of_vertices();


        for (int j = 0; j < Points.size() ; j++)
        {
            vertex_map[Points[j]] = j;
        }
        //! EdgeWeightChanger:

        boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, G);


        Triangulation::Finite_edges_iterator iter;
        for(iter =  T.finite_edges_begin();
                iter != T.finite_edges_end();
                iter++)
        {
            Triangulation::Triangulation_data_structure::Edge e = *iter;
            Triangulation::Triangulation_data_structure::Cell_handle c = e.first;
            int i = e.second;
            int j = e.third;
            Point pa = c->vertex(i)->point();
            Point pb = c->vertex(j)->point();
            int p = vertex_map[pa];
            int q = vertex_map[pb];

            //  Graph::add_edge(Tmp,vertex_map[pa],vertex_map[pb]);
            //     G[e].distance = d;



            //! Handle edges here:
            double newDistance = sqrt(CGAL::squared_distance(pa,pb));
            //auto additionStat = Graph::add
            int pn = componentMap[p];
            int qn = componentMap[q];
            if (pn != qn)
            {
                auto additionstat = boost::add_edge(pn,qn,G);
                auto ee = additionstat.first;
                if (additionstat.second)
                {
                    weightmap[ee] = newDistance;
                    G[ee].distance = newDistance;
                }
                else
                {
                    if (newDistance < weightmap[ee])
                    {
                        weightmap[ee] = newDistance;
                        G[ee].distance = newDistance;
                    }

                }

            }
        }

    }

    void FormSimplexTree(MyGraphType & G, Gudhi::Simplex_tree<> & complexx)
    {

        Filtration_value zero = 0;
        //! Insert vertices:
        auto a = boost::vertices(G);
        for (auto it = a.first; it != a.second; it++)
        {
            int v = *it;
            typeVectorVertex abc = {v};
            complexx.insert_simplex_and_subfaces(abc, 0);
        }

        //! Insert edges:
        auto b = boost::edges(G);
        for (auto it = b.first; it != b.second; it++)
        {
            int v = boost::source(*it,G);
            int w = boost::target(*it,G);
            typeVectorVertex edgee = {v,w};
            double ddd = G[*it].distance;
            // Filtration_value tmp = G[*it].distance;
            complexx.insert_simplex_and_subfaces(edgee,ddd);
        }

        //! Complete the thing:
         complexx.expansion(2);
        std::cout << "Vertex count: " << complexx.num_vertices() << std::endl;
        std::cout << "Simplices: " << complexx.num_simplices() << std::endl;


        std::cout << "**************************************************************" << std::endl;
        std::cout << "strict graph G { " << std::endl;
        for (auto f_simplex : complexx.filtration_simplex_range())
        {
            std::cout << "   " << "[" << complexx.filtration(f_simplex) << "] ";
            for (auto vertex : complexx.simplex_vertex_range(f_simplex))
            {
                std::cout << static_cast<int>(vertex) << " -- ";
            }
            std::cout << ";" << std::endl;
        }
        std::cout << "}" << std::endl;
        std::cout << "**************************************************************" << std::endl;



    }





    void Run(std::vector<Point> & points, double fEpsilon, double clusterEpsilon, MyGraphType & G, std::string colorout, Gudhi::Simplex_tree<> & complexx)
    {
        std::vector<size_t> clusterCorrespondance(points.size());
        std::vector<size_t> componentMap(points.size());
        std::vector<size_t> leaders;
        std::vector<size_t> resultingMap;
        int num = 0;
        ComputeClustering(points,clusterCorrespondance,leaders,fEpsilon,clusterEpsilon,num,componentMap,colorout);
        DefineMapToMajorElement(num, componentMap, clusterCorrespondance, resultingMap);
        ComputeDeluanayTriangulation( G, clusterCorrespondance,resultingMap, componentMap, points,num);
        FormSimplexTree(G,complexx);
    }

protected:

private:

};


//  void ClusterAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2, std::unordered_map<int, std::vector<int>> & clusters)
//    {
//        //! 1) We create tree
//        std::cout << "Treestage" << std::endl;
//        std::vector<size_t> oldfromnew;
//        Tree* kTree = new Tree(in, oldfromnew);
//        //! 2a) We compute the KDE
//        std::vector<double> f(in.n_cols);
//        std::vector<int> vn(in.n_cols);
//        DensityComputator<MetricType, MatType, DensityType, TreeType> calcf(kTree);
//        calcf.ComputeDensity(f, vn, epsilon);
//       // arma::mat dg = in.t();
//        //! 2b) Rescaling
//        Rescale(f, t);
//        //! 3) We compute new MST
//        EngineCluster<MetricType, MatType, TreeType> finalCalc(kTree,oldfromnew, false);
//        finalCalc.ComputeAMST(out,f,epsilon2,clusters);
//        kTree = NULL;
//        delete kTree;
//
//    }
#endif // CLUSTERPERSISTENCE_H

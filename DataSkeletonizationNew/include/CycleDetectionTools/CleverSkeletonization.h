#ifndef CLEVERSKELETONIZATION_H
#define CLEVERSKELETONIZATION_H

#include "AmstStat.h"
#include <mlpack/methods/emst/edge_pair.hpp>
#include <mlpack/methods/emst/union_find.hpp>
#include <mlpack/methods/emst/dtb.hpp>
#include <mlpack/methods/emst/dtb_rules.hpp>
#include "EpsilonClusterMST/EpsilonClusteringRules.h"
#include "unordered_map"
#include "GeneralConvertor.h"
#include "Graph.h"

//! The idea is to build mst on clusters
//! In first step we map every element to the largest value in the epsilon neighborhood its located
//! In the second we compute MSt on existing clusterization with one exception
template<
    typename MetricType,
    typename MatType,
    template<typename TreeMetricType,
             typename TreeStatType,
             typename TreeMatType> class TreeType>
class CleverSkeletonization
{
public:
    typedef TreeType<MetricType, AmstStat, MatType> Tree;
private:
    //! The function f ()
    //  std::vector<double>* f;
    //! Permutations of points during tree building.
    std::vector<size_t> oldFromNew;
    //! Pointer to the root of the tree.
    Tree* tree;
    //! Reference to the data (this is what should be used for accessing data).
    const MatType& data;
    //! Indicates whether or not we "own" the tree.
    bool ownTree;
    //! IndicatesEpsilonClusterMST/ whether or not O(n^2) naive mode will be used.
    bool naive;
    //! List of edge nodes.
    arma::Col<size_t> neighborsInComponent;
    //! List of edge nodes.
    arma::Col<size_t> neighborsOutComponent;
    //! List of edge distances.
    arma::vec neighborsDistances;
    //! Total distance of the tree.
    double totalDist;
    //! The instantiated metric.
    MetricType metric;
    //! Specialmode
    bool specialmode;
    //! New function required for this method.
    std::vector<size_t> clusterCorrespondance;


    struct SortEdgesHelper
    {
        bool operator()(const mlpack::emst::EdgePair& pairA, const mlpack::emst::EdgePair& pairB)
        {
            return (pairA.Distance() < pairB.Distance());
        }
    } SortFun;

public:
    CleverSkeletonization(  Tree* tree, std::vector<size_t> & ofn, bool naive = false,
                            const MetricType metric = MetricType()
                         ) :
        tree(tree),
        data(tree->Dataset()),
        ownTree(false),
        naive(naive),
        totalDist(0.0),
        metric(metric),
        specialmode(true),
        oldFromNew(ofn)

    {
        clusterCorrespondance.resize(data.n_cols);
        neighborsInComponent.set_size(data.n_cols);
        neighborsOutComponent.set_size(data.n_cols);
        neighborsDistances.set_size(data.n_cols);
        neighborsDistances.fill(DBL_MAX);

    }

    //! Works currently only for trees which contain all the explicit points in leafs and all non-leaf elements do not contain any points
    //! Because we were too lazy to fix this xD


    void RecursiveComputator(std::vector<double> & f, Tree* node)
    {
        if (node->IsLeaf())
        {
            for (int i = 0; i < node->NumPoints(); i++)
            {
                double fl = f[node->Point(i)];
                node->Stat().maxF() = std::max(node->Stat().maxF(),fl);
                if (fl < node->Stat().minF())
                {
                    node->Stat().setMinFIndex(node->Point(i));
                }
                node->Stat().minF() = std::min(node->Stat().minF(),fl);
            }
        }
        else
        {
            for (int i = 0; i < node->NumChildren(); i++)
            {
                Tree* child = &(node->Child(i));
                RecursiveComputator(f, child);
                node->Stat().maxF() = std::max(node->Stat().maxF(),child->Stat().maxF());
                if (child->Stat().minF() < node->Stat().minF())
                {
                    node->Stat().setMinFIndex(child->Stat().returnMinFIndex());
                }
                node->Stat().minF() = std::min(node->Stat().minF(),child->Stat().minF());
            }
        }
    }




    //! We updated the Graph from MyGraphType to AbstractGraphType
    //! (because we dont need the extra information provided by the MyGraphType parameters

    int CompleteClusterizaiton(int & num, std::vector<size_t> & componentMap)
    {
        //! Add all the elements to the graph
        AbstractGraphType G;
        for (int i = 0; i < data.n_cols; i++)
        {
            boost::add_vertex(G);
        }
        int componentsFormed = 0;

        //! Identify two types of points: Global roots and local roots  AND
        //! Add all connections to the graph

        std::vector<size_t> endPointDetector(data.n_cols);
        std::vector<size_t> localRoot(data.n_cols);
        std::vector<size_t> localRootSize(data.n_cols);
        std::list<size_t> listOfRoots;
        for (int i = 0; i < clusterCorrespondance.size(); i++)
        {
            int j = clusterCorrespondance[i];
            boost::add_edge(i,j,G);
            localRoot[j] = 1;
            localRootSize[j]++;
            if(i == j)
            {
                endPointDetector[j] = 1;
                listOfRoots.push_back(i);
            }
            else
            {
                componentsFormed++;
            }
        }
        //! Calculate the map for connected components
        // Commented out, because we want to own this information:
        // std::vector<int> componentMap(num_vertices(G));
        componentMap.resize(num_vertices(G));
        num = boost::connected_components(G, componentMap.data());
        std::cout << "Currently, componentMap size: " << componentMap.size() << std::endl;
        std::cout << "Currenlty num size: " << num << std::endl;
        std::vector<size_t> componentsToRootsMap(num);
        for (auto it = listOfRoots.begin(); it != listOfRoots.end(); it++)
        {
            componentsToRootsMap[componentMap[*it]] = *it;
        }

        //! Update the cluster correspondance map and print some info

        std::vector<size_t> rootComponentSize(data.n_cols);
        std::vector<std::vector<size_t>> subClusterSizes(data.n_cols);

        for (int i = 0; i < clusterCorrespondance.size(); i++)
        {
            clusterCorrespondance[i] = componentsToRootsMap[componentMap[i]];
        }

//        //! Do stupid calculations and useless printing here:
//        for (int i = 0; i < clusterCorrespondance.size(); i++)
//        {
//            rootComponentSize[clusterCorrespondance[i]]++;
//            if (localRoot[i] == 1)
//            {
//                subClusterSizes[clusterCorrespondance[i]].push_back(localRootSize[i]);
//            }
//        }
//        for (auto it = listOfRoots.begin(); it != listOfRoots.end(); it++)
//        {
//            std::cout << "Size of the cluster: " << rootComponentSize[*it] << " | Split into: ";
//            for (int j = 0; j < subClusterSizes[*it].size(); j++)
//            {
//                std::cout << subClusterSizes[*it][j] << ", ";
//            }
//            std::cout << std::endl;
//
//        }
//        std::cout << "THe number" << num << std::endl;
        return componentsFormed;

    }



    void ComputeClusterization(std::vector<size_t> & ClusterCorrespondanceFinal, std::vector<double> & f, double epsilon, int & num, std::vector<size_t> & componentMap)
    {
        mlpack::Timer::Start("amst/amst_computation");
        totalDist = 0; // Reset distance.

        //! Rearrange f if required.
        std::vector<double>* fp = &f;
        if (ownTree && mlpack::tree::TreeTraits<Tree>::RearrangesDataset)
        {
            fp = new std::vector<double>(data.n_cols);
            for (int i = 0; i < data.n_cols; i++)
            {
                const size_t nindex = oldFromNew[i];
                (*fp)[nindex] = f[i];
            }
        }
        //! Setup the tree for the computations
        RecursiveComputator(*fp,tree);
//        typedef RulesRangeAMST<MetricType, Tree> RuleType2;

        typedef EpsilonClusteringRules<MetricType, Tree> RuleType;
        RuleType rules(data, neighborsDistances, metric, *fp, clusterCorrespondance, epsilon);

        //! Perform the clustering and remember to connect the clusters in the unionfind data-structure

        typename Tree::template DualTreeTraverser<RuleType> traverser(rules);
        traverser.Traverse(*tree, *tree);

        std::vector<size_t> componentMapTmp;
        int componentsFormed = CompleteClusterizaiton(num,componentMapTmp);

        //! Allocated components:
        ClusterCorrespondanceFinal.resize(clusterCorrespondance.size());
        for (int i = 0; i < this-> clusterCorrespondance.size(); i++ )
        {
            int nindex = oldFromNew[i];
            ClusterCorrespondanceFinal[nindex] = oldFromNew[clusterCorrespondance[i]];
            componentMap[nindex] = componentMapTmp[i];
        }


        fp = NULL;
        delete fp;
    }



};

#endif // CLEVERSKELETONIZATION

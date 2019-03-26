#ifndef DUALTREEAMST_H
#define DUALTREEAMST_H

#include "AmstStat.h"
#include "RulesAMST.h"
#include "RulesRangeAMST.h"
#include <mlpack/methods/emst/edge_pair.hpp>
#include <mlpack/methods/emst/union_find.hpp>
#include "RulesRAMSTM.h"

//! Call the tree constructor that does mapping.
template<typename TreeType, typename MatType>
TreeType* BuildTree(
    MatType&& dataset,
    std::vector<size_t>& oldFromNew,
    const typename std::enable_if<
    mlpack::tree::TreeTraits<TreeType>::RearrangesDataset>::type* = 0)
{
    return new TreeType(std::forward<MatType>(dataset), oldFromNew);
}

//! Call the tree constructor that does not do mapping.
template<typename TreeType, typename MatType>
TreeType* BuildTree(
    MatType&& dataset,
    const std::vector<size_t>& /* oldFromNew */,
    const typename std::enable_if<
    !mlpack::tree::TreeTraits<TreeType>::RearrangesDataset>::type* = 0)
{
    return new TreeType(std::forward<MatType>(dataset));
}

template<
    typename MetricType,
    typename MatType,
    template<typename TreeMetricType,
             typename TreeStatType,
             typename TreeMatType> class TreeType>
class DualTreeAMST
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

    //! Indicates whether or not O(n^2) naive mode will be used.
    bool naive;

    //! Edges.
    std::vector<mlpack::emst::EdgePair> edges; // We must use vector with non-numerical types.

    //! Connections.
    mlpack::emst::UnionFind connections;

    //! List of edge nodes.
    arma::Col<size_t> neighborsInComponent;
    //! List of edge nodes.
    arma::Col<size_t> neighborsOutComponent;
    //! List of edge distances.
    arma::vec neighborsDistances;
    //! Tiebreak distances
    arma::vec tiebreakDistances;

    //! Total distance of the tree.
    double totalDist;
    //! The instantiated metric.
    MetricType metric;
//! Specialmode
    bool specialmode;

    struct SortEdgesHelper
    {
        bool operator()(const mlpack::emst::EdgePair& pairA, const mlpack::emst::EdgePair& pairB)
        {
            return (pairA.Distance() < pairB.Distance());
        }
    } SortFun;

public:



    DualTreeAMST( const MatType& dataset,
                  const bool naive = false,
                  const MetricType metric = MetricType()) :
        tree(naive ? NULL : BuildTree<Tree>(dataset, oldFromNew)),
        data(naive ? dataset : tree->Dataset()),
        ownTree(!naive),
        naive(naive),
        connections(dataset.n_cols),
        totalDist(0.0),
        metric(metric),
        specialmode(false)
    {
        edges.reserve(data.n_cols - 1); // Set size.

        neighborsInComponent.set_size(data.n_cols);
        neighborsOutComponent.set_size(data.n_cols);
        neighborsDistances.set_size(data.n_cols);
        neighborsDistances.fill(DBL_MAX);
        tiebreakDistances.set_size(data.n_cols);
        tiebreakDistances.fill(DBL_MAX);
    }

    DualTreeAMST(  Tree* tree,
                   const MetricType metric = MetricType(),
                   bool naive = false) :
        tree(tree),
        data(tree->Dataset()),
        ownTree(false),
        naive(naive),
        connections(data.n_cols),
        totalDist(0.0),
        metric(metric),
        specialmode(false)
    {
        edges.reserve(data.n_cols - 1); // Fill with EdgePairs.

        neighborsInComponent.set_size(data.n_cols);
        neighborsOutComponent.set_size(data.n_cols);
        neighborsDistances.set_size(data.n_cols);
        neighborsDistances.fill(DBL_MAX);
        tiebreakDistances.set_size(data.n_cols);
        tiebreakDistances.fill(DBL_MAX);
    }
    DualTreeAMST(  Tree* tree, std::vector<size_t> & ofn, bool naive = false,
                   const MetricType metric = MetricType()
                ) :
        tree(tree),
        data(tree->Dataset()),
        ownTree(false),
        naive(naive),
        connections(data.n_cols),
        totalDist(0.0),
        metric(metric),
        specialmode(true),
        oldFromNew(ofn)

    {
        edges.reserve(data.n_cols - 1); // Fill with EdgePairs.

        neighborsInComponent.set_size(data.n_cols);
        neighborsOutComponent.set_size(data.n_cols);
        neighborsDistances.set_size(data.n_cols);
        neighborsDistances.fill(DBL_MAX);
        tiebreakDistances.set_size(data.n_cols);
        tiebreakDistances.fill(DBL_MAX);
    }

    ~DualTreeAMST()
    {
        if (ownTree)
            delete tree;
    }
    void AddEdge(
        const size_t e1,
        const size_t e2,
        const double distance)
    {
        mlpack::Log::Assert((distance >= 0.0),
                            "DualTreeAMST::AddEdge(): distance cannot be negative.");

        if (e1 < e2)
            edges.push_back(mlpack::emst::EdgePair(e1, e2, distance));
        else
            edges.push_back(mlpack::emst::EdgePair(e2, e1, distance));
    }
    void AddAllEdges()
    {
        for (size_t i = 0; i < data.n_cols; i++)
        {
            size_t component = connections.Find(i);
            size_t inEdge = neighborsInComponent[component];
            size_t outEdge = neighborsOutComponent[component];
            if (connections.Find(inEdge) != connections.Find(outEdge))
            {
                // totalDist = totalDist + dist;
                // changed to make this agree with the cover tree code
                totalDist += neighborsDistances[component];
                AddEdge(inEdge, outEdge, neighborsDistances[component]);
                connections.Union(inEdge, outEdge);
            }
        }
    }


    //! Just a blank method currently
    void ComputeFInTree()
    {
        //   std::vector<double>* pf = &f;
    }
    void RecursiveComputator(std::vector<double> & f, Tree* node)
    {
        if (node->IsLeaf())
        {
            for (int i = 0; i < node->NumPoints(); i++)
            {
                double fl = f[node->Point(i)];
                node->Stat().maxF() = std::max(node->Stat().maxF(),fl);
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
                node->Stat().minF() = std::min(node->Stat().minF(),child->Stat().minF());
            }
        }
    }


    void EmitResults(arma::mat& results)
    {
        // Sort the edges.
        std::sort(edges.begin(), edges.end(), SortFun);

        mlpack::Log::Assert(edges.size() == data.n_cols - 1);
        results.set_size(3, edges.size());

        // Need to unpermute the point labels.
        if ((!naive && ownTree && mlpack::tree::TreeTraits<Tree>::RearrangesDataset)||(specialmode == true))
        {
            for (size_t i = 0; i < (data.n_cols - 1); i++)
            {
                // Make sure the edge list stores the smaller index first to
                // make checking correctness easier
                size_t ind1 = oldFromNew[edges[i].Lesser()];
                size_t ind2 = oldFromNew[edges[i].Greater()];

                if (ind1 < ind2)
                {
                    edges[i].Lesser() = ind1;
                    edges[i].Greater() = ind2;
                }
                else
                {
                    edges[i].Lesser() = ind2;
                    edges[i].Greater() = ind1;
                }

                results(0, i) = edges[i].Lesser();
                results(1, i) = edges[i].Greater();
                results(2, i) = edges[i].Distance();
            }
        }
        else
        {
            for (size_t i = 0; i < edges.size(); i++)
            {
                results(0, i) = edges[i].Lesser();
                results(1, i) = edges[i].Greater();
                results(2, i) = edges[i].Distance();
            }
        }
    }

    void CleanupHelper(Tree* tree)
    {
        // Reset the statistic information.
        tree->Stat().MaxNeighborDistance() = DBL_MAX;
        tree->Stat().MinNeighborDistance() = DBL_MAX;
        tree->Stat().Bound() = DBL_MAX;

        // Recurse into all children.
        for (size_t i = 0; i < tree->NumChildren(); ++i)
            CleanupHelper(&tree->Child(i));

        // Get the component of the first child or point.  Then we will check to see
        // if all other components of children and points are the same.
        const int component = (tree->NumChildren() != 0) ?
                              tree->Child(0).Stat().ComponentMembership() :
                              connections.Find(tree->Point(0));

        // Check components of children.
        for (size_t i = 0; i < tree->NumChildren(); ++i)
            if (tree->Child(i).Stat().ComponentMembership() != component)
                return;

        // Check components of points.
        for (size_t i = 0; i < tree->NumPoints(); ++i)
            if (connections.Find(tree->Point(i)) != size_t(component))
                return;

        // If we made it this far, all components are the same.
        tree->Stat().ComponentMembership() = component;
    }

    void Cleanup()
    {
        for (size_t i = 0; i < data.n_cols; i++)
        {
            neighborsDistances[i] = DBL_MAX;
            tiebreakDistances[i] = DBL_MAX;
            }

        if (!naive)
            CleanupHelper(tree);
    }

    void ComputeAMST(arma::mat& results, std::vector<double> & f)
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


        typedef RulesAMST<MetricType, Tree> RuleType;
        RuleType rules(data, connections, neighborsDistances, neighborsInComponent,
                       neighborsOutComponent, metric, *fp);
        while (edges.size() < (data.n_cols - 1))
        {
            if (naive)
            {
                // Full O(N^2) traversal.
                for (size_t i = 0; i < data.n_cols; ++i)
                    for (size_t j = 0; j < data.n_cols; ++j)
                        rules.BaseCase(i, j);
            }
            else
            {
                typename Tree::template DualTreeTraverser<RuleType> traverser(rules);
                traverser.Traverse(*tree, *tree);
            }

            AddAllEdges();

            Cleanup();

            mlpack::Log::Info << edges.size() << " edges found so far." << std::endl;
            if (!naive)
            {
                //! We are too lazy to implement this currently
                //   mlpack::Log::Info << rules.BaseCases() << " cumulative base cases." << std::endl;
                //  mlpack::Log::Info << rules.Scores() << " cumulative node combinations scored."
                //            << std::endl;
            }
        }

        mlpack::Timer::Stop("amst/amst_computation");

        EmitResults(results);

        mlpack::Log::Info << "Total (augemented) spanning tree length: " << totalDist << std::endl;
        fp = NULL;
        delete fp;
    }
    void ComputeRangeAMST(arma::mat& results, std::vector<double> & f, double epsilon)
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

        typedef RulesRangeAMST<MetricType, Tree> RuleType2;
        RuleType2 rules(data, connections, neighborsDistances, tiebreakDistances, neighborsInComponent,
                       neighborsOutComponent, metric, *fp, epsilon);
        while (edges.size() < (data.n_cols - 1))
        {
            if (naive)
            {
                // Full O(N^2) traversal.
                for (size_t i = 0; i < data.n_cols; ++i)
                    for (size_t j = 0; j < data.n_cols; ++j)
                        rules.BaseCase(i, j);
            }
            else
            {
                typename Tree::template DualTreeTraverser<RuleType2> traverser(rules);
                traverser.Traverse(*tree, *tree);
            }

            AddAllEdges();

            Cleanup();

            mlpack::Log::Info << edges.size() << " edges found so far." << std::endl;
            if (!naive)
            {
                //! We are too lazy to implement this currently
                //   mlpack::Log::Info << rules.BaseCases() << " cumulative base cases." << std::endl;
                //  mlpack::Log::Info << rules.Scores() << " cumulative node combinations scored."
                //            << std::endl;
            }
        }

        mlpack::Timer::Stop("amst/amst_computation");

        EmitResults(results);

        mlpack::Log::Info << "Total (augemented) spanning tree length: " << totalDist << std::endl;
        fp = NULL;
        delete fp;
    }
     void ComputeRangeAMSTModifiedTwo(arma::mat& results, std::vector<double> & f, double epsilon)
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

        typedef RulesRAMSTM<MetricType, Tree> RuleType3;
        RuleType3 rules(data, connections, neighborsDistances, neighborsInComponent,
                       neighborsOutComponent, metric, *fp, epsilon);
        while (edges.size() < (data.n_cols - 1))
        {
            if (naive)
            {
                // Full O(N^2) traversal.
                for (size_t i = 0; i < data.n_cols; ++i)
                    for (size_t j = 0; j < data.n_cols; ++j)
                        rules.BaseCase(i, j);
            }
            else
            {
                typename Tree::template DualTreeTraverser<RuleType3> traverser(rules);
                traverser.Traverse(*tree, *tree);
            }

            AddAllEdges();

            Cleanup();

            mlpack::Log::Info << edges.size() << " edges found so far." << std::endl;
            if (!naive)
            {
                //! We are too lazy to implement this currently
                //   mlpack::Log::Info << rules.BaseCases() << " cumulative base cases." << std::endl;
                //  mlpack::Log::Info << rules.Scores() << " cumulative node combinations scored."
                //            << std::endl;
            }
        }

        mlpack::Timer::Stop("amst/amst_computation");

        EmitResults(results);

        mlpack::Log::Info << "Total (augemented) spanning tree length: " << totalDist << std::endl;
        fp = NULL;
        delete fp;
    }




protected:

};

#endif // DUALTREEAMST_H

#ifndef EPSILONCLUSTERINGRULES_H
#define EPSILONCLUSTERINGRULES_H
#include <mlpack/methods/emst/union_find.hpp>
#include <mlpack/prereqs.hpp>
#include <mlpack/core/tree/traversal_info.hpp>

template<typename MetricType, typename TreeType>
class EpsilonClusteringRules
{
public:
    EpsilonClusteringRules(const arma::mat& dataSet,
                           arma::vec& neighborsDistances,
                           MetricType& metric,
                           const std::vector<double> & f,
                           std::vector<size_t> & correspondance,
                           double epsilon
                          ): dataSet(dataSet),
        neighborsDistances(neighborsDistances),
        metric(metric),
        baseCases(0),
        scores(0),
        f(f),
        correspondance(correspondance),
        epsilon(epsilon)
    {}
    double BaseCase(const size_t queryIndex,
                    const size_t referenceIndex)
    {
        // Check if the points are in the same component at this iteration.
        // If not, return the distance between them.  Also, store a better result as
        // the current neighbor, if necessary.

        double newUpperBound = 0;

        ++baseCases;
        double distance = metric.Evaluate(dataSet.col(queryIndex),
                                          dataSet.col(referenceIndex));

        // double thevalue = (std::min(f[queryIndex],f[referenceIndex]))*distance;

        if (distance < epsilon)
        {
            double thevalue = f[referenceIndex];
            if (thevalue < neighborsDistances[queryIndex])
            {
                neighborsDistances[queryIndex] = thevalue;
                correspondance[queryIndex] = referenceIndex;
            }
        }




        return newUpperBound;
    }

    double Score(const size_t queryIndex,TreeType& referenceNode)
    {


        const arma::vec queryPoint = dataSet.unsafe_col(queryIndex);
        const double distance = referenceNode.MinDistance(queryPoint);
        if (distance > epsilon)
        {
            return DBL_MAX;
        }

        // If all the points in the reference node are farther than the candidate
        // nearest neighbor for the query's component, we prune.
        return neighborsDistances[queryIndex] < referenceNode.Stat().minF()
               ? DBL_MAX : referenceNode.Stat().minF() ;
    }

    double Rescore(const size_t queryIndex,
                   TreeType& /* referenceNode */,
                   const double oldScore)
    {
        // We don't need to check component membership again, because it can't
        // change inside a single iteration.
        return (oldScore > neighborsDistances[queryIndex])
               ? DBL_MAX : oldScore;
    }
    double Score(TreeType& queryNode, TreeType& referenceNode)
    {
        // If all the queries belong to the same component as all the references
        // then we prune.

        //! In case max distance is under epsilon, we can go directly to the lowest level.
        //! We are however too lazy to implement it

        const double distance = queryNode.MinDistance(referenceNode);

        //const double valuev = queryNode.Stat().minF()/2 + referenceNode.Stat().minF()/2;

        //! Case when we can do rapid descend:

        int indxOfMinNode = referenceNode.Stat().returnMinFIndex();
        const double maxDistanceToMinNode = queryNode.MaxDistance(dataSet.col(indxOfMinNode));

        if (maxDistanceToMinNode < epsilon)
        {
            //! For every point in queryNode:

            for (size_t t = 0; t < queryNode.NumDescendants(); t++)
            {
                size_t indexOfDescendant = queryNode.Descendant(t);
                BaseCase(indexOfDescendant, indxOfMinNode);
            }
            return DBL_MAX;
        }

        //! Case when the distance is too big

        if (distance > epsilon)
        {
            return DBL_MAX;
        }

        //! Other cases


        double fvalue = referenceNode.Stat().minF();
        // const double bound = CalculateBound(queryNode);
        // If all the points in the reference node are farther than the candidate
        // nearest neighbor for all queries in the node, we prune.
        ++scores;

        return distance;
    }

    double Rescore(TreeType& queryNode,
                   TreeType& /* referenceNode */,
                   const double oldScore) const
    {
        return oldScore;
    }

    typedef typename mlpack::tree::TraversalInfo<TreeType> TraversalInfoType;

    const TraversalInfoType& TraversalInfo() const
    {
        return traversalInfo;
    }
    TraversalInfoType& TraversalInfo()
    {
        return traversalInfo;
    }

protected:

private:

    TraversalInfoType traversalInfo;

    //! The function f
    const std::vector<double> & f;

    //! The data points.
    const arma::mat& dataSet;

    //! The distance to the candidate nearest neighbor for each component.
    arma::vec& neighborsDistances;

    //! The instantiated metric.
    MetricType& metric;

    //! The number of base cases calculated.
    size_t baseCases;
    //! The number of node combinations that have been scored.
    size_t scores;
    //! Epsilon
    double epsilon;
    //!  The new function
    std::vector<size_t>& correspondance;


};

#endif // RULESRANGEAMST_H

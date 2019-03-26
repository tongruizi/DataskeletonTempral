#ifndef RULESRANGEAMST_H
#define RULESRANGEAMST_H
#include <mlpack/methods/emst/union_find.hpp>
#include <mlpack/prereqs.hpp>

template<typename MetricType, typename TreeType>
class RulesRangeAMST
{
public:
    RulesRangeAMST(const arma::mat& dataSet,
                   mlpack::emst::UnionFind& connections,
                   arma::vec& neighborsDistances,
                   arma::vec& tiebreakDistances,
                   arma::Col<size_t>& neighborsInComponent,
                   arma::Col<size_t>& neighborsOutComponent,
                   MetricType& metric,
                   const std::vector<double> & f,
                   double epsilon
                  ): dataSet(dataSet),
        connections(connections),
        neighborsDistances(neighborsDistances),
        tiebreakDistances(tiebreakDistances),
        neighborsInComponent(neighborsInComponent),
        neighborsOutComponent(neighborsOutComponent),
        metric(metric),
        baseCases(0),
        scores(0),
        f(f),
        epsilon(epsilon)
    {}
    double BaseCase(const size_t queryIndex,
                    const size_t referenceIndex)
    {
        // Check if the points are in the same component at this iteration.
        // If not, return the distance between them.  Also, store a better result as
        // the current neighbor, if necessary.
        double newUpperBound = -1.0;

        // Find the index of the component the query is in.
        size_t queryComponentIndex = connections.Find(queryIndex);
        size_t referenceComponentIndex = connections.Find(referenceIndex);
        if (queryComponentIndex != referenceComponentIndex)
        {
            ++baseCases;
            double distance = metric.Evaluate(dataSet.col(queryIndex),
                                              dataSet.col(referenceIndex));

            // double thevalue = (std::min(f[queryIndex],f[referenceIndex]))*distance;

            if (distance < epsilon)
            {
                double thevalue = f[referenceIndex];
                double tiebreakvalue = f[queryIndex];
                // f[queryIndex]/2+
                mlpack::Log::Assert(queryIndex != referenceIndex);
                if (neighborsOutComponent[queryComponentIndex] == referenceIndex)
                {
                    if (tiebreakvalue < tiebreakDistances[queryIndex])
                    {
                    tiebreakDistances[queryIndex] = tiebreakvalue;
                    neighborsDistances[queryComponentIndex] = thevalue;
                    //! Just to make sure :D
                    neighborsInComponent[queryComponentIndex] = queryIndex;
                    neighborsOutComponent[queryComponentIndex] = referenceIndex;
                    }

                }
                else if (thevalue < neighborsDistances[queryComponentIndex])
                {
                    neighborsDistances[queryComponentIndex] = thevalue;
                    neighborsInComponent[queryComponentIndex] = queryIndex;
                    neighborsOutComponent[queryComponentIndex] = referenceIndex;
                }
            }
        }

        if (newUpperBound < neighborsDistances[queryComponentIndex])
            newUpperBound = neighborsDistances[queryComponentIndex];

        mlpack::Log::Assert(newUpperBound >= 0.0);

        return newUpperBound;
    }

    double Score(const size_t queryIndex,TreeType& referenceNode)
    {
        size_t queryComponentIndex = connections.Find(queryIndex);

        // If the query belongs to the same component as all of the references,
        // then prune.  The cast is to stop a warning about comparing unsigned to
        // signed values.
        if (queryComponentIndex ==
                (size_t) referenceNode.Stat().ComponentMembership())
            return DBL_MAX;

        const arma::vec queryPoint = dataSet.unsafe_col(queryIndex);
        const double distance = referenceNode.MinDistance(queryPoint);
        if (distance > epsilon)
        {
            return DBL_MAX;
        }

        // If all the points in the reference node are farther than the candidate
        // nearest neighbor for the query's component, we prune.
        return neighborsDistances[queryComponentIndex] < referenceNode.Stat().minF()
               ? DBL_MAX : referenceNode.Stat().minF() ;
    }

    double Rescore(const size_t queryIndex,
                   TreeType& /* referenceNode */,
                   const double oldScore)
    {
        // We don't need to check component membership again, because it can't
        // change inside a single iteration.
        return (oldScore > neighborsDistances[connections.Find(queryIndex)])
               ? DBL_MAX : oldScore;
    }
    double Score(TreeType& queryNode, TreeType& referenceNode)
    {
        // If all the queries belong to the same component as all the references
        // then we prune.
        if ((queryNode.Stat().ComponentMembership() >= 0) &&
                (queryNode.Stat().ComponentMembership() ==
                 referenceNode.Stat().ComponentMembership()))
            return DBL_MAX;

        const double distance = queryNode.MinDistance(referenceNode);
        //const double valuev = queryNode.Stat().minF()/2 + referenceNode.Stat().minF()/2;

        if (distance > epsilon)
        {

            return DBL_MAX;
        }
        double fvalue = referenceNode.Stat().minF();

        const double bound = CalculateBound(queryNode);

        // If all the points in the reference node are farther than the candidate
        // nearest neighbor for all queries in the node, we prune.
        ++scores;

        return (bound < fvalue) ? DBL_MAX : distance;
    }

    double Rescore(TreeType& queryNode,
                   TreeType& /* referenceNode */,
                   const double oldScore) const
    {
        const double bound = CalculateBound(queryNode);
        return (oldScore > bound) ? DBL_MAX : oldScore;
    }
    double CalculateBound(TreeType& queryNode) const
    {
        double worstPointBound = -DBL_MAX;
        double worstChildBound = -DBL_MAX;

        // Now, find the best and worst point bounds.
        for (size_t i = 0; i < queryNode.NumPoints(); ++i)
        {
            const size_t pointComponent = connections.Find(queryNode.Point(i));
            const double bound = neighborsDistances[pointComponent];

            if (bound > worstPointBound)
                worstPointBound = bound;
        }
        // Find the best and worst child bounds.
        for (size_t i = 0; i < queryNode.NumChildren(); ++i)
        {
            const double maxBound = queryNode.Child(i).Stat().MaxNeighborDistance();
            if (maxBound > worstChildBound)
                worstChildBound = maxBound;

        }
        // Now calculate the actual bounds.
        const double worstBound = std::max(worstPointBound, worstChildBound);

        // Update the relevant quantities in the node.
        queryNode.Stat().MaxNeighborDistance() = worstBound;
        return queryNode.Stat().Bound();
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

    //! Stores the tree structure so far
    mlpack::emst::UnionFind& connections;

    //! The distance to the candidate nearest neighbor for each component.
    arma::vec& neighborsDistances;

    //! Closest f.
    arma::vec& tiebreakDistances;


    //! The index of the point in the component that is an endpoint of the
    //! candidate edge.
    arma::Col<size_t>& neighborsInComponent;

    //! The index of the point outside of the component that is an endpoint
    //! of the candidate edge.
    arma::Col<size_t>& neighborsOutComponent;

    // arma::Col<size_t>&


    //! The instantiated metric.
    MetricType& metric;

    //! The number of base cases calculated.
    size_t baseCases;
    //! The number of node combinations that have been scored.
    size_t scores;
    //! Epsilon
    double epsilon;


};

#endif // RULESRANGEAMST_H

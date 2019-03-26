#ifndef DENSITYESTIMATIONRULES_H
#define DENSITYESTIMATIONRULES_H

#include <vector>
#include <mlpack/prereqs.hpp>

template<typename MetricType,
         typename MatType,
         typename FunctionType,
         typename TreeType>

class DensityEstimationRules
{
public:
// WRITE THIS FUNCTION Page x in the paper.
    DensityEstimationRules(MatType & referenceSetd, MatType & querySetd, std::vector<double> & one, std::vector<int> & visits, double ep):
        referenceSet(referenceSetd),
        querySet(querySetd),
        results(one),
        visitNumber(visits),
        epsilon(ep),
        baseCases(0),
        scores(0)
    {

    }
    double BaseCase(const size_t queryIndex, const size_t referenceIndex)
    {
        if ((lastQueryIndex == queryIndex) && (lastReferenceIndex == referenceIndex))
        {
            return 0;
        }
        double distt = MetricType::Evaluate(querySet.unsafe_col(queryIndex), referenceSet.unsafe_col(referenceIndex));
        double fvalue = FunctionType::Evaluate(distt);
        baseCases++;
        lastQueryIndex = queryIndex;
        lastReferenceIndex = referenceIndex;
        results[queryIndex] = results[queryIndex] + fvalue;
        visitNumber[queryIndex] = visitNumber[queryIndex] + 1;
        return fvalue;

    }
    double Score(const size_t queryIndex, TreeType& referenceNode)
    {
        double ep = epsilon;
        double mxvalue = referenceNode.MaxDistance((querySet).col(queryIndex));
        double minvalue = referenceNode.MinDistance((querySet).col(queryIndex));
        double wvalue = FunctionType::Evaluate((minvalue)) - FunctionType::Evaluate((mxvalue));
        //  referenceSet.n_cols
        if (wvalue <= ep)
        {
            int multiplier = referenceNode.NumDescendants();
            arma::colvec first = ((querySet).col(queryIndex));
            arma::colvec second;
            referenceNode.Center(second);
            results[queryIndex] = results[queryIndex] + multiplier * FunctionType::Evaluate(MetricType::Evaluate(first,second));
            visitNumber[queryIndex] = visitNumber[queryIndex] + multiplier;
            return DBL_MAX;

        }
        return minvalue;
    }
     double Rescore(const size_t /*index */, TreeType& /* referenceNode */, const double oldScore)
    {
        return oldScore;
    }


    // IMPLEMENT THE REQUIRED METHOD FOR SCORE
    double Score(TreeType & queryNode, TreeType & referenceNode)
    {
        double ep = epsilon;
        double mxvalue = referenceNode.MaxDistance(queryNode);
        double minvalue = referenceNode.MinDistance(queryNode);
        double wvalue = FunctionType::Evaluate((minvalue)) - FunctionType::Evaluate(mxvalue);
        if (wvalue <= ep)
        {
            int num_descendants = queryNode.NumDescendants();
            for (int i = 0; i < num_descendants; i++)
            {
                int indx = queryNode.Descendant(i);
                int multiplier = referenceNode.NumDescendants();
                arma::colvec first = querySet.col(indx);
                arma::colvec second;
                referenceNode.Center(second);
                results[indx] = results[indx] + multiplier * FunctionType::Evaluate(MetricType::Evaluate(first,second));
                visitNumber[indx] = visitNumber[indx] + multiplier;

            }
            return DBL_MAX;
        }
        return minvalue;

    }

    double Rescore(TreeType& /* queryNode */, TreeType& /* referenceNode */, const double oldScore)
    {
        return oldScore;
    }

    void TestMethodAdditionToResults(int k, int n)
    {
        results[k] = n;

    }
    double returnbaseCases()
    {
        return baseCases;
    }
    double returnScores()
    {
        return scores;
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
    double epsilon;
    std::vector<double>& results;
    std::vector<int>& visitNumber;
    MatType& referenceSet;
    MatType& querySet;
    size_t lastQueryIndex;
    size_t lastReferenceIndex;
    size_t baseCases;
    size_t scores;
    //TraversalInfoType traversalInfo;

    TraversalInfoType traversalInfo;

    //std::vector<double> sums;
};

#endif // DENSITYESTIMATIONRULES_H

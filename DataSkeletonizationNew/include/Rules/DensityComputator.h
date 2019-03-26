#ifndef DENSITYCOMPUTATOR_H
#define DENSITYCOMPUTATOR_H

// This class is used to start computation of density estimation.
#include <mlpack/prereqs.hpp>
#include <mlpack/core.hpp>
#include <mlpack/core/metrics/lmetric.hpp>
#include <mlpack/core/tree/binary_space_tree.hpp>
#include "DensityEstimationRules.h"
#include "ExponentialDensity.h"
#include "KDEStats.h"
#include "TreeCreator.h"
#include "AmstStat.h"

template<typename MetricType,
         typename MatType,
         typename FunctionType,
         template<typename TreeMetricType,
                  typename TreeStatType,
                  typename TreeMatType> class TreeType>
class DensityComputator
{
    typedef TreeType<MetricType,AmstStat, MatType> Tree;
    // typedef
    //typedef TreeType Tree;


private:
    //! Mappings to old reference indices (used when this object builds trees).
    std::vector<size_t> oldFromNewReferences;
    //! Our tree:
    Tree* referenceTree;
    //! Our possibly rearranged referenceset:
    MatType* referenceSet;
    size_t baseCases;
    size_t scores;
    //! Owns the tree
    bool treeOwner;

public:

    DensityComputator(MatType referenceSet):
        referenceTree(TreeCreator::BuildTree<Tree>(std::move(referenceSet), oldFromNewReferences)),
        referenceSet(new MatType(referenceTree->Dataset())),
        baseCases(0),
        scores(0),
        treeOwner(true)
    {

    }
    DensityComputator(Tree* referenceTree):
        referenceTree(referenceTree),
        referenceSet(&referenceTree->Dataset()),
        treeOwner(false),
        baseCases(0),
        scores(0)
    {
    }
    ~DensityComputator()
    {
        if (treeOwner && referenceTree)
            delete referenceTree;

    }

    void ComputeDensity(std::vector<double> & values, std::vector<int> & visitNumber, double epsilon)
    {
        mlpack::Timer::Start("KDE/Computation");
        std::vector<double>* valuesPtr = & values;
        std::vector<int>* visitNumberPtr = & visitNumber;

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            valuesPtr = new std::vector<double>(referenceSet->n_cols);
            visitNumberPtr = new std::vector<int>(referenceSet->n_cols);
        }

        typedef DensityEstimationRules<MetricType,MatType,FunctionType,Tree> RuleType;
        RuleType rules(*referenceSet, *referenceSet, *valuesPtr,*visitNumberPtr, epsilon);

        typename Tree::template DualTreeTraverser<RuleType> traverser(rules);
        traverser.Traverse(*referenceTree, *referenceTree);
        //     baseCases = rules.returnBaseCases();
        //    scores = rules.returnScores();

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            values.clear();
            values.resize(referenceSet->n_cols);
            visitNumber.clear();
            visitNumber.resize(referenceSet->n_cols);
            for (int i = 0; i < values.size(); i++)
            {
                const size_t nindex = oldFromNewReferences[i];
                values[nindex] = (*valuesPtr)[i];
                visitNumber[nindex] = (*visitNumberPtr)[i];
            }
        }
        //  else
        //   {
        //      for (int i = 0; i < values.size(); i++)
        //      {
        //          values[i] = (*valuesPtr)[i];
        //         visitNumber[i] = (*visitNumberPtr)[i];
        //      }

        //    }
        mlpack::Timer::Stop("KDE/Computation");
        delete valuesPtr;
        valuesPtr = NULL;
        delete visitNumberPtr;
        visitNumberPtr = NULL;

    }
    //! We have disabled timer for now (to avoid extra struggle)

    void ComputeDensity(std::vector<double> & values, std::vector<int> & visitNumber, double epsilon, Tree* queryTree)
    {
        //    mlpack::Timer::Start("KDE/Computation");
        const MatType& querySet = queryTree->Dataset();
        std::vector<double>* valuesPtr = & values;
        std::vector<int>* visitNumberPtr = & visitNumber;

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            valuesPtr = new std::vector<double>(querySet->n_cols);
            visitNumberPtr = new std::vector<int>(querySet->n_cols);
        }

        typedef DensityEstimationRules<MetricType,MatType,FunctionType,Tree> RuleType;
        RuleType rules(*referenceSet, queryTree->Dataset(), *valuesPtr,*visitNumberPtr, epsilon);

        typename Tree::template DualTreeTraverser<RuleType> traverser(rules);
        traverser.Traverse(*queryTree, *referenceTree);
        //     baseCases = rules.returnBaseCases();
        //    scores = rules.returnScores();

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            values.clear();
            values.resize(valuesPtr->size());
            visitNumber.clear();
            visitNumber.resize(visitNumberPtr->size());
            for (int i = 0; i < values.size(); i++)
            {
                const size_t nindex = oldFromNewReferences[i];
                values[nindex] = (*valuesPtr)[i];
                visitNumber[nindex] = (*visitNumberPtr)[i];
            }
        }
        //  else
        //   {
        //      for (int i = 0; i < values.size(); i++)
        //      {
        //          values[i] = (*valuesPtr)[i];
        //         visitNumber[i] = (*visitNumberPtr)[i];
        //      }

        //    }
        //  mlpack::Timer::Stop("KDE/Computation");
        delete valuesPtr;
        valuesPtr = NULL;
        delete visitNumberPtr;
        visitNumberPtr = NULL;
    }

    void ComputeDensity(std::vector<double> & values, std::vector<int> & visitNumber, double epsilon, const MatType& querySet)
    {
        //    mlpack::Timer::Start("KDE/Computation");
        std::vector<double>* valuesPtr = & values;
        std::vector<int>* visitNumberPtr = & visitNumber;

        //! This is important vector
        std::vector<size_t> oldFromNewQueries;

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            valuesPtr = new std::vector<double>(querySet.n_cols);
            visitNumberPtr = new std::vector<int>(querySet.n_cols);
        }

        typedef DensityEstimationRules<MetricType,MatType,FunctionType,Tree> RuleType;

        Tree* queryTree = TreeCreator::BuildTree<Tree>(querySet, oldFromNewQueries);

        RuleType rules(*referenceSet, queryTree->Dataset(), *valuesPtr,*visitNumberPtr, epsilon);


        typename Tree::template DualTreeTraverser<RuleType> traverser(rules);
        traverser.Traverse(*queryTree, *referenceTree);
        //     baseCases = rules.returnBaseCases();
        //    scores = rules.returnScores();

        if(mlpack::tree::TreeTraits<Tree>::RearrangesDataset && treeOwner)
        {
            values.clear();
            values.resize(valuesPtr->size());
            visitNumber.clear();
            visitNumber.resize(visitNumberPtr->size());
            for (int i = 0; i < values.size(); i++)
            {
                const size_t nindex = oldFromNewQueries[i];
                values[nindex] = (*valuesPtr)[i];
                visitNumber[nindex] = (*visitNumberPtr)[i];
            }
        }
        //  else
        //   {
        //      for (int i = 0; i < values.size(); i++)
        //      {
        //          values[i] = (*valuesPtr)[i];
        //         visitNumber[i] = (*visitNumberPtr)[i];
        //      }

        //    }
        //  mlpack::Timer::Stop("KDE/Computation");
        delete valuesPtr;
        valuesPtr = NULL;
        delete visitNumberPtr;
        visitNumberPtr = NULL;

    }




    double ReturnZero()
    {
        return 0;
    }

protected:



};

#endif // COMPUTEDENSITYESTIMATION_H

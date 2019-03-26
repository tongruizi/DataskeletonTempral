#ifndef AMSTCOMPUTATOR_H
#define AMSTCOMPUTATOR_H

#include "DualTreeAMST.h"
#include "DensityComputator.h"
#include "EpsilonClusterMST/EngineClusterMST.h"

template<
    typename MetricType,
    typename MatType,
    typename DensityType,
    template<typename TreeMetricType,
             typename TreeStatType,
             typename TreeMatType> class TreeType>
class AMSTComputator
{
public:
    typedef TreeType<MetricType, AmstStat, MatType> Tree;

    AMSTComputator() {}

    void DebugKDE(const std::vector<double> & oldKD, std::vector<size_t> & tr, arma::mat & data, std::string output)
    {
        std::vector<double> newKD(oldKD.size());
        for (int i = 0; i < newKD.size(); i++)
        {
            newKD[tr[i]] = oldKD[i];
        }
        GeneralConvertor::MatInfoToFile(output, data, newKD);

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
    void PerformAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t)
    {
        std::cout << "The first element: " << in(0,0) << std::endl;
        std::cout << "Succesfully launched" << std::endl;
        //! 1) We create tree
        std::vector<size_t> oldfromnew;
        Tree* kTree = new Tree(in, oldfromnew);
        // Tree* kTree = TreeCreator::BuildTree<Tree>(std::move(in), oldfromnew);
        std::cout << "Tree allocated succefully " << std::endl;
        //! 2a) We compute the KDE
        std::vector<double> f(in.n_cols);
        std::vector<int> vn(in.n_cols);
        DensityComputator<MetricType, MatType, DensityType, TreeType> calcf(kTree);
        std::cout << "Before the computation" << std::endl;
        calcf.ComputeDensity(f, vn, epsilon);
        std::cout << "Visit numbers: " << vn[0] << std::endl;
        std::cout << "Computed KDE succefully " << std::endl;
        arma::mat dg = in.t();
        //! 2b) Rescaling
        DebugKDE(f, oldfromnew, dg, "/home/yury/Dropbox/MlPackTraining/KernelDensityEstimation/outputs/amst/beforenormdebug.csv");
        Rescale(f, t);
        DebugKDE(f, oldfromnew, dg, "/home/yury/Dropbox/MlPackTraining/KernelDensityEstimation/outputs/amst/afternormdebug.csv");

        //! 3) We compute the tree
        DualTreeAMST<MetricType, MatType, TreeType> finalCalc(kTree,oldfromnew, true);
        finalCalc.ComputeAMST(out, f);
        std::cout << "Computed the aMST succefully " << std::endl;
        // delete kTree;
    }
    void PerformRangeAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2)
    {
        //! 1) We create tree
        std::cout << "Treestage" << std::endl;
        std::vector<size_t> oldfromnew;
        Tree* kTree = new Tree(in, oldfromnew);
        // Tree* kTree = TreeCreator::BuildTree<Tree>(std::move(in), oldfromnew);
        //! 2a) We compute the KDE
        std::cout << "KDEstage" << std::endl;
        std::vector<double> f(in.n_cols);
        std::vector<int> vn(in.n_cols);
        DensityComputator<MetricType, MatType, DensityType, TreeType> calcf(kTree);
        calcf.ComputeDensity(f, vn, epsilon);
        arma::mat dg = in.t();
        //! 2b) Rescaling
        Rescale(f, t);
        //! 3) We compute the tree
        std::cout << "AMSTstage" << std::endl;
        DualTreeAMST<MetricType, MatType, TreeType> finalCalc(kTree,oldfromnew, false);
        //finalCalc.ComputeRangeAMST(out, f, epsilon2);
        finalCalc.ComputeRangeAMSTModifiedTwo(out,f,epsilon2);
        // delete kTree;
        kTree = NULL;
        delete kTree;
        std::cout << "Completed" << std::endl;
    }
    void ClusterAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2, std::unordered_map<int, std::vector<int>> & clusters)
    {
        //! 1) We create tree
        std::cout << "Treestage" << std::endl;
        std::vector<size_t> oldfromnew;
        Tree* kTree = new Tree(in, oldfromnew);
        //! 2a) We compute the KDE
        std::vector<double> f(in.n_cols);
        std::vector<int> vn(in.n_cols);
        DensityComputator<MetricType, MatType, DensityType, TreeType> calcf(kTree);
        calcf.ComputeDensity(f, vn, epsilon);
       // arma::mat dg = in.t();
        //! 2b) Rescaling
        Rescale(f, t);
        //! 3) We compute new MST
        EngineCluster<MetricType, MatType, TreeType> finalCalc(kTree,oldfromnew, false);
        finalCalc.ComputeAMST(out,f,epsilon2,clusters);
        kTree = NULL;
        delete kTree;

    }

    virtual ~AMSTComputator() {}

protected:

private:
};

#endif // AMSTCOMPUTATOR_H

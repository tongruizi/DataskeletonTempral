#ifndef COMPUTATOR_H
#define COMPUTATOR_H

#include <mlpack/core.hpp>
#include "ExponentialDensity.h"
#include "DensityComputator.h"
#include "GeneralConvertor.h"
#include "ToSptTransformer.h"
#include "DiMoSCLauncher.h"

//#include "DiMoSC.cpp"
//! This guy requires some fixing: #include "DiMoSC.cpp"

// Row,Column

//   k(0,0) = 0;
//     k(1,0) = 2;
class MorseReconstructionComputator
{
public:
    typedef  DensityComputator<mlpack::metric::EuclideanDistance, arma::mat, ExponentialDensity, mlpack::tree::KDTree> ourDensityComputator;
    MorseReconstructionComputator() {}

    // out,data,proportion,densityE,parC
    static void ComputeWithProportions(MyGraphType & G, arma::mat & data, int msd, double epsilon, double parameter)
    {

        arma::vec MaxValues(data.n_rows);
        arma::vec MinValues(data.n_rows);
        MinValues.fill(std::numeric_limits<double>::max());

        arma::vec Difference(data.n_rows);
        for (int i = 0; i < data.n_cols; i++)
        {
            for (int j = 0; j < data.n_rows; j++)
            {
                MaxValues(j) = std::max(data(j,i), MaxValues(j));
                MinValues(j) = std::min(data(j,i), MinValues(j));
            }
        }

        double smallestDifference = std::numeric_limits<double>::max();
        for (int j = 0; j < data.n_rows; j++)
        {
            Difference(j) = MaxValues(j)-MinValues(j);
            smallestDifference = std::min(smallestDifference,Difference(j));
        }

        std::vector<int> cellsizes(data.n_rows);

        for(int j = 0; j < data.n_rows; j++)
        {
            double ratioo = (double)msd * ((double) Difference(j) / (double) smallestDifference);
            cellsizes[j] = ceil(ratioo);
        }

        CoreComputation(G,cellsizes,epsilon,data,MaxValues,MinValues,parameter);

    }

    static void Compute(MyGraphType & G, arma::mat & data, std::vector<int> & cellsizes, double epsilon, double parameter)
    {


        arma::vec MaxValues(data.n_rows);
        arma::vec MinValues(data.n_rows);
        MinValues.fill(std::numeric_limits<double>::max());
        for (int i = 0; i < data.n_cols; i++)
        {
            for (int j = 0; j < data.n_rows; j++)
            {
                MaxValues(j) = std::max(data(j,i), MaxValues(j));
                MinValues(j) = std::min(data(j,i), MinValues(j));
            }
        }
        CoreComputation(G,cellsizes,epsilon,data,MaxValues,MinValues,parameter);

    }

    // /home/yury/LocalTests/DebugDebug
    // static void PrintMatrixToFile(std::string & w, arma::mat & matrix)

    static void TransformDataAndLaunchDimosc(MyGraphType & G, arma::mat & data, arma::Cube<int> & IndexationStater, std::vector<double> & f, int x, int y, int z, double parameter)
    {
        arma::Mat<int> tetra;
        ToSptTransformer::buildTetraGrid(x,y,z,tetra,IndexationStater);
        // std::cout << tetra << std::endl;
        arma::Mat<int> tri;
        ToSptTransformer::buildTriFromTetra(tetra,tri);

        arma::Mat<int> edge;
        ToSptTransformer::buildEdgeFromTri(tri,edge);

        DiMoSCLauncher::Launch(G,parameter,data,f,edge,tri);
        // MyGraphType & G, double ve_delta, arma::mat & data, std::vector<double> & f, arma::Mat<int> & edges, arma::Mat<int> & triangles
    }



    static void CoreComputation(MyGraphType & G, std::vector<int> & cellsizes, double epsilon, arma::mat & data, arma::vec & MaxValues, arma::vec & MinValues, double parameter)
    {

        int querySetSize = 1;

        arma::Cube<int> IndexationStater(cellsizes[0],cellsizes[1],cellsizes[2]);

        for (int k = 0; k < cellsizes.size(); k++)
        {
            querySetSize = querySetSize * cellsizes[k];
        }
        arma::mat querySet(data.n_rows, querySetSize);
        //  std::cout << querySet << std::endl;
        std::vector<int> counter(cellsizes.size(),0);
        for (int i = 0; i < querySetSize; i++)
        {
            IndexationStater(counter[0],counter[1],counter[2]) = i;
            for (int j = 0; j <  cellsizes.size() ; j++)
            {
                double fraction1 = (double) counter[j]/ (double) (cellsizes[j] - 1);
                double ddd = MinValues(j)*(1-fraction1) + MaxValues(j)*(fraction1);
                querySet(j,i) = ddd;
            }
            int tmp = cellsizes.size() - 1;
            while(counter[tmp] == cellsizes[tmp] - 1)
            {
                counter[tmp] = 0;
                tmp--;
            }
            if (i != querySetSize-1)
            {
                counter[tmp] = counter[tmp] + 1;
            }
        }

        //! Use the algorithms:

        std::vector<double> SumVector(querySet.n_cols);
        std::vector<int> visitNumber(querySet.n_cols);
        ourDensityComputator DC(data);
        //   ComputeDensity(std::vector<double> & values, std::vector<int> & visitNumber, double epsilon, Tree* queryTree)

        DC.ComputeDensity(SumVector, visitNumber, epsilon,querySet);


        //! Convert Back:
        //! if needed   querySet = querySet.t();

        // GeneralConvertor::normalizedOut(pathout,cellsizes[0],cellsizes[1],cellsizes[2],SumVector);
        //   GeneralConvertor::ScalarToFile(pathout,cellsizes[0],cellsizes[1],cellsizes[2],SumVector);
        // GeneralConvertor::MatInfoToFile(pathout,querySet,SumVector);
        TransformDataAndLaunchDimosc(G, querySet, IndexationStater,SumVector,cellsizes[0],cellsizes[1],cellsizes[2],parameter);
    }


protected:

private:
};

#endif // COMPUTATOR_H

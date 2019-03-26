#ifndef TESTPERFORMER_H
#define TESTPERFORMER_H

#include "AMSTComputator.h"
#include <mlpack/core.hpp>
#include "CloudGenerator.h"
#include "Computation.h"
#include "BranchDetection.h"
#include "GraphGeneration.h"


class TestPerformer
{
public:
    TestPerformer() {}

    static void RunTestsQuick(std::string & outFolder, int gn, int n, double scale, double minangle, int k, double epc, double epd, double epmst)
    {

        double t = 2;
        double timesumm = 0;
        AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, ExponentialDensity,mlpack::tree::KDTree> comp;
        for (int i = 0; i < k; i++)
        {
            std::unordered_map<int, std::vector<int>> clusterz;
            std::string numberI =  std::to_string(i);
            std::string outputName = "outputAMST" + numberI + ".vtk";
            std::string outputPath = outFolder + outputName;

            arma::mat cloud;
            arma::mat outputGraph;

            //! Use this for standart testing:
              CloudGenerator::generatePointsStandartInArmaMat(gn, n,minangle,scale,epc,cloud);
            //! Use this for stupid triangle testing:
//            MyGraphType A;
//            std::list<Point> pointList;
//            GraphGeneration::TriangleGraph(A, scale);
//            CloudGenerator::generatePoints(n,A,epc,pointList);
//            GeneralConvertor::ListToMat(pointList, cloud);

            //! Continue from here
            arma::mat trcloud = cloud.t();

            // ClusterAMSTComputation(arma::mat & in, arma::mat & out, double epsilon, double t, double epsilon2)
            comp.ClusterAMSTComputation(trcloud,outputGraph,epd,t,epmst,clusterz);
            GeneralConvertor::MSTToVTK(trcloud,outputGraph, outputPath);

            //! Print out the clusters:
            //! We wont print clusters for a while, mostly because we are too lazy to fix dimensionality issues

//            int numfirst = 0;
//            std::cout << "Just before cluster number thing" << std::endl;
//            for (auto it = clusterz.begin(); it != clusterz.end(); it++)
//            {
//            std::cout << "Number of element in: " << it->first << " is: " << (it->second).size() << std::endl;
//
//            }
//            for (auto it = clusterz.begin(); it != clusterz.end(); it++)
//            {
//                std::cout << "Here we are " << std::endl;
//                std::string debugpath = outFolder + "debug" + std::to_string(numfirst) + ".csv";
//                std::vector<double> simpleFunction(cloud.n_rows);
//                int elemnumber = (it->second).size();
//                std::cout << "Before theCloud" << std::endl;
//                arma::mat theCloud(elemnumber,3);
//               // theCloud = theCloud.()
//                for (int j = 0; j < elemnumber; j++)
//                {
//                    for (int p = 0; p < 3; p++)
//                    {
//                    theCloud(j,p) = cloud(j,p);
//                    }
//                }
//                std::cout << "After theCloud" << std::endl;
//                simpleFunction[it->first] = 1;
//                //arma::mat newcloud =  theCloud.t();
//                std::cout << theCloud << std::endl;
//                GeneralConvertor::MatInfoToFile(debugpath, theCloud, simpleFunction);
//                numfirst++;
//            }

        }
    }

    static void RunTests(std::string & outFolder, int gn, int n, double scale, double minangle, int k, double epc, double epd, double epmst)
    {
        double t = 2;
        double timesumm = 0;
        //! The method:
        AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, ExponentialDensity,mlpack::tree::KDTree> comp;
        for (int i = 0; i < k; i++)
        {
            std::string numberI =  std::to_string(i);
            std::string outputName = "output" + numberI + ".vtk";
            std::string outputPath = outFolder + outputName;
            arma::mat cloud;
            arma::mat outputGraph;
            CloudGenerator::generatePointsStandartInArmaMat(gn, n,minangle,scale,epc,cloud);
            arma::mat trcloud = cloud.t();
            mlpack::Timer::Start("test/Computation1" + numberI);
            comp.PerformRangeAMSTComputation(trcloud,outputGraph,epd,t,epmst);
            mlpack::Timer::Stop("test/Computation1" + numberI);

            GeneralConvertor::MSTToVTK(trcloud,outputGraph, outputPath);
            auto gtm = mlpack::Timer::Get("test/Computation1" + numberI);
            long tmtm = gtm.count();
            timesumm = timesumm + tmtm;
            MyGraphType G;
            MyGraphType K;
            GeneralConvertor::ArmaMatToGraph(G,outputGraph,trcloud);
            BranchDetection::SimplifyIt(G,K,50.0,"","pure");

            std::string ww = outFolder + "finaloutput" + numberI + ".vtk";
            GeneralConvertor::GraphToVtk(ww,K);

            std::cout << "Time for test " << i << " : " << tmtm << std::endl;

        }
        std::cout << "On AVERAGE: " << timesumm/k << std::endl;
    }
    static void RunCyclicTests(std::string & outFolder, int n, double scale, int k, double epc, double epd, double epmst)
    {
        double t = 2;
        double timesumm = 0;
        AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, ExponentialDensity,mlpack::tree::KDTree> comp;

        for (int i = 0; i < k; i++)
        {
            arma::mat cloud;
            MyGraphType A;
            std::list<Point> pointList;
            std::string numberI =  std::to_string(i);
            std::string outputName = "output" + numberI + ".vtk";
            std::string outputPath = outFolder + outputName;
            std::string outputPath2 = outFolder + "FinalOutput" + numberI + ".vtk";
            GraphGeneration::TriangleGraph(A, scale);
            CloudGenerator::generatePoints(n,A,epc,pointList);
            GeneralConvertor::ListToMat(pointList, cloud);

            arma::mat outputGraph;
            arma::mat trcloud = cloud.t();
            mlpack::Timer::Start("test/Computation1" + numberI);
            comp.PerformRangeAMSTComputation(trcloud,outputGraph,epd,t,epmst);
            mlpack::Timer::Stop("test/Computation1" + numberI);

            GeneralConvertor::MSTToVTK(trcloud,outputGraph, outputPath);
            auto gtm = mlpack::Timer::Get("test/Computation1" + numberI);
            long tmtm = gtm.count();
            timesumm = timesumm + tmtm;
            MyGraphType G;
            MyGraphType K;
            GeneralConvertor::ArmaMatToGraph(G,outputGraph,trcloud);
            BranchDetection::SimplifyIt(G,K,50.0,"","pure");
            GeneralConvertor::GraphToVtk(outputPath2,K);

            std::cout << "Time for test " << i << " : " << tmtm << std::endl;

        }
        std::cout << "On AVERAGE: " << timesumm/k << std::endl;

    }

protected:

private:
};

#endif // TESTPERFORMER_H

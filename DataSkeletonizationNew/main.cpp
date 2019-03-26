double EPS_compare = 0.0000001;

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "ExponentialDensity.h"
#include "LinearDensity.h"
#include "FunctionEvaluator.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/methods/emst/dtb.hpp>
#include "DensityComputator.h"
#include "GeneralConvertor.h"
#include "NaiveMethods.h"
#include "DensityEstimationRules.h"
#include "TreeExperements.h"
#include <mlpack/core/tree/statistic.hpp>
#include "AMSTComputator.h"
#include <ctime>
#include "TestPerformer.h"
#include <cstdlib>
#include "GradientDescendTester.h"
//! Mapper and alphareeb includes:
#include "AlphaReeb_Launcher.h"
#include "CorrectEndTypeMeasure.h"
#include "ExplicitMeasurer.h"
#include "AbstractMeasurer.h"
#include "Mapper_Launcher.h"

#include "SingleStar.h"
#include "ClassicDistanceMeasure.h"
#include "controller.h"
#include "PrintToFile.h"
#include "NumberOfVertexMeasure.h"
#include "CorrectTypeMeasure.h"
#include "AsKLauncher.h"
#include "DoubleStar.h"

#include <boost/filesystem.hpp>
#include "RealCloudCollection.h"
#include "FastMSTLauncher.h"
#include "DualTreeDistanceMeasure.h"
#include "MeanSquareDistanceMeasure.h"
#include "DumbAlgorithm.h"
#include "CloudTypePrinterAlgorithm.h"
#include "SophisticatedPrinter.h"
#include "AmstLauncher.h"
#include "SegmentDistance.h"

#include "MorseLauncher.h"
#include "MetricGraphRec/MetricRec.hpp"




// Convenience.
using namespace mlpack;
using namespace mlpack::neighbor;
using namespace mlpack::metric;
// Convenience.

using namespace std;






void ControllerTest()
{
    std::string qq = "/home/yury/LocalTests/Test1/data.csv";
    std::string folder = "/home/yury/LocalTests/Test1/Outputs/";
    double mappercluster = 1.75; // 1.75
    double alpha = 20;
//! We initilize a star:
    SingleStar star1(M_PI/3,3,1500,5,100,5,"Star3");
    SingleStar star2(M_PI/3,4,2000,5,100,5,"Star4");
    SingleStar star8(M_PI/6,8,4000,5,100,1,"Star8");
    DoubleStar dstar(M_PI/3,4,4,4000,5,100,10,"DoubleStar");

    // RealCloudCollection theRealClouds("Real","/home/yury/Dropbox/UnileverData/XYZ_Files/");
    //theRealClouds.SetCorrectnessOfGraphs("/home/yury/Dropbox/UnileverData/titles.txt");
//RealCloudCollection theRealClouds("Real","/home/yury/Dropbox/UnileverData/XYZ_Files/");
//SingleStar star2(M_PI/3,4,100,1500,5,100,10,"Star4");
//! We initilize filewriter:
    PrintToFile printer(folder);
//! We initilize mapper:
    Mapper_Parameters param(15, 0.5, "Distance", mappercluster,mappercluster);
    Mapper_Launcher thelaunch(param);
    AlphaReeb_Parameters AlphaParam(alpha, 1);
    AlphaReeb_Launcher alphaLaunch(AlphaParam,10);
    thelaunch.setTimePrecision(2);
    alphaLaunch.setTimePrecision(2);
    std::string setting = "sd";
    //  AsKLauncher AskAlgorithm(2.0,1.3,1.3,setting,"AsK");
    AsKLauncher AskAlgorithm(20.0,1.3,1.3,"MSTAVG","AsKMsTAvG");
    AskAlgorithm.setTimePrecision(3);

    //  AskAlgorithm.addPostRunner(&printer);
//  thelaunch.addPostRunner(&printer);
    //  alphaLaunch.addPostRunner(&printer);

    //! Add AMST:

    AmstLauncher lwp(0.01,1,4,15.0);
    lwp.setTimePrecision(3);
    lwp.addPostRunner(&printer);

//! Add MST

    FastMSTLauncher mstComputator;
    mstComputator.setTimePrecision(3);
    mstComputator.addPostRunner(&printer);

//! We initilize controller:
    controller control(qq);
//! We add algorithm, star to controller
//   void addAlgorithm(AbstractAlgorithm* k);
//   void addCloud(generatable* k);
    //control.addAlgorithm(&thelaunch);
    // control.addAlgorithm(&alphaLaunch);
    // control.addAlgorithm(&AskAlgorithm);
    control.addAlgorithm(&mstComputator);
// control.addAlgorithm(&lwp);
//control.addCloud(&star1);
//control.addCloud(&star2);
    control.addCloud(&star8);
//control.addCloud(&theRealClouds);
//control.addCloud(&dstar);
//control.addCloud(&star2);
//control.addCloud(&theRealClouds);
//! Initilize distance error measure:
    ClassicDistanceMeasure distanceMeasure(2);
    CorrectEndTypeMeasure endTypeMeasure(2);
    NumberOfVertexMeasure vertexMeasure(2);
    CorrectTypeMeasure TypeMeasure(2);
    DualTreeDistanceMeasure newDistanceMeasure(2);
    MeanSquareDistanceMeasure ndm(2);

    control.addMeasurer(distanceMeasure);
    control.addMeasurer(endTypeMeasure);
    control.addMeasurer(vertexMeasure);
    control.addMeasurer(TypeMeasure);
    control.addMeasurer(newDistanceMeasure);
    // control.addMeasurer(ndm);
    // control.addMeasurer(newDistanceMeasure);
    control.BeginTestRun();

}

void ControllerTestRealDataSimple()
{
//! Defining main algorithm parameters:
    double mappercluster = 3.0; // 1.75
    double alpha = 20;

//! Defining all required file paths:
    std::string CSVFile = "/home/yury/Dropbox/UnileverData/FirstTestOfInterface/data2.csv";
    std::string Outputfolder = "/home/yury/Dropbox/UnileverData/FirstTestOfInterface/Outputs2/";
    std::string Inputfolder = "/home/yury/Dropbox/UnileverData/FirstTestOfInterface/SecondClusters/";
    std::string ManualGraphInfoFile = "/home/yury/Dropbox/UnileverData/FirstTestOfInterface/StarInfo2.txt";


//! Initilizing RealCloudCollection
    RealCloudCollection theRealClouds("Real",Inputfolder);
    std::cout << "Initilized RealCloudCollection" << std::endl;
    theRealClouds.SetCorrectnessOfGraphs(ManualGraphInfoFile);
    std::cout << "Succesfully imputed file " << std::endl;

    std::cout << "Initilization of RealCloudCollection succeful" << std::endl;

//! Initilizing printer:
    PrintToFile printer(Outputfolder);


//! Initilizing Mapper Algorithm
    Mapper_Parameters param(15, 0.5, "Distance", mappercluster,mappercluster);
    Mapper_Launcher MapperAlgorithm(param);
    MapperAlgorithm.setTimePrecision(2);
    MapperAlgorithm.addPostRunner(&printer);

//! Initilizing AsK Algorithm
    std::string setting = "pure";
    AsKLauncher AskAlgorithm(50.0,1.05,1.5,setting,"AsK");
    AskAlgorithm.setTimePrecision(2);
    AskAlgorithm.addPostRunner(&printer);

//! Initilize MST
    FastMSTLauncher mstComputator;
    mstComputator.setTimePrecision(2);
    mstComputator.addPostRunner(&printer);

//! Initilizing the controller:

    controller control(CSVFile);
    control.addCloud(&theRealClouds);
//    control.addAlgorithm(&MapperAlgorithm);
    control.addAlgorithm(&mstComputator);
    control.addAlgorithm(&AskAlgorithm);


//! Initilizing measurers:

    ClassicDistanceMeasure distanceMeasure(2);
    CorrectEndTypeMeasure endTypeMeasure(2);
    NumberOfVertexMeasure vertexMeasure(2);
    CorrectTypeMeasure TypeMeasure(2);

//! Adding measurers

    control.addMeasurer(distanceMeasure);
    control.addMeasurer(endTypeMeasure);
    control.addMeasurer(vertexMeasure);
    control.addMeasurer(TypeMeasure);

//! Run the boy
    control.BeginTestRun();

}

void DistanceBetweenLineSegments()
{
    Segment l(Point(0,0,10), Point(1,0,10));
    Segment q(Point(0,1,0), Point(0,1,0));
    std::cout << CGAL::squared_distance(l, q) << std::endl;

}

void PrintingStarCloudsTest()
{
    DumbAlgorithm pwee;
    SophisticatedPrinter p("/home/yury/LocalTests/TestOnSynthetical/");
    pwee.addPostRunner(&p);
    SingleStar star3(M_PI/3,3,1500,5,100,100,"Star3");
    SingleStar star4(M_PI/3,4,2000,5,100,100,"Star4");
    SingleStar star5(M_PI/3,5,2500,5,100,100,"Star5");
    SingleStar star6(M_PI/6,6,3000,5,100,100,"Star6");
    SingleStar star7(M_PI/3,7,3500,5,100,100,"Star7");
    SingleStar star8(M_PI/3,8,4000,5,100,100,"Star8");
    std::vector<CloudTypePrinterAlgorithm*> q;
    CloudTypePrinterAlgorithm n3("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star3Info.txt");
    CloudTypePrinterAlgorithm n4("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star4Info.txt");
    CloudTypePrinterAlgorithm n5("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star5Info.txt");
    CloudTypePrinterAlgorithm n6("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star6Info.txt");
    CloudTypePrinterAlgorithm n7("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star7Info.txt");
    CloudTypePrinterAlgorithm n8("/home/yury/LocalTests/TestOnSynthetical/TheBosses/Star8Info.txt");
    q.push_back(&n3);
    q.push_back(&n4);
    q.push_back(&n5);
    q.push_back(&n6);
    q.push_back(&n7);
    q.push_back(&n8);
    controller control("/home/yury/LocalTests/TestOnSynthetical/dummyfile.txt");
    control.addAlgorithm(&pwee);
    control.addCloud(&star3);
    control.addCloud(&star4);
    control.addCloud(&star5);
    control.addCloud(&star6);
    control.addCloud(&star7);
    control.addCloud(&star8);
    control.BeginTestRun();
    control.FlushClouds(q);

//#include "CloudTypePrinterAlgorithm.h"



}
//

void ControllerTest1337()
{
  //  std::string qq = "/home/liudi/Dropbox/LocalTests/TwoAlgorithmTest/data.csv";
  std::string qq = "/home/liudi/Dropbox/LocalTests/OutputMetricGraphRec/debugcolorActual.csv";
    std::string folder = "/home/liudi/Dropbox/LocalTests/TwoAlgorithmTest/Outputs/";
//! We initilize a star:
    SingleStar star8(M_PI/6,8,4000,5,100,1,"Star8");


//! We initilize filewriter:
    PrintToFile printer(folder);

// int proportion, double densityE, double parC
    MorseLauncher morseAlgorithm(10,0.01, 0.10);
    morseAlgorithm.setTimePrecision(3);
    morseAlgorithm.addPostRunner(&printer);



//! Add MST
//! Initilize MST
    FastMSTLauncher mstComputator;
    mstComputator.setTimePrecision(2);
    mstComputator.addPostRunner(&printer);

//! We initilize controller:
    controller control(qq);
//! We add algorithm, star to controller
//   void addAlgorithm(AbstractAlgorithm* k);
//   void addCloud(generatable* k);
    //control.addAlgorithm(&thelaunch);
    // control.addAlgorithm(&alphaLaunch);
    // control.addAlgorithm(&AskAlgorithm);
    control.addAlgorithm(&morseAlgorithm);
    control.addAlgorithm(&mstComputator);
// control.addAlgorithm(&lwp);
//control.addCloud(&star1);
//control.addCloud(&star2);
    control.addCloud(&star8);
//control.addCloud(&theRealClouds);
//control.addCloud(&dstar);
//control.addCloud(&star2);
//control.addCloud(&theRealClouds);
//! Initilize distance error measure:
    ClassicDistanceMeasure distanceMeasure(2);
    CorrectEndTypeMeasure endTypeMeasure(2);
    NumberOfVertexMeasure vertexMeasure(2);
    CorrectTypeMeasure TypeMeasure(2);
    DualTreeDistanceMeasure newDistanceMeasure(2);
    MeanSquareDistanceMeasure ndm(2);

    control.addMeasurer(distanceMeasure);
    control.addMeasurer(endTypeMeasure);
    control.addMeasurer(vertexMeasure);
    control.addMeasurer(TypeMeasure);
    // control.addMeasurer(ndm);
    // control.addMeasurer(newDistanceMeasure);
    control.BeginTestRun();




}


void MetriccontrollerTest()
{
  //  std::string qq = "/home/liudi/Dropbox/LocalTests/TwoAlgorithmTest/data.csv";
  std::string qq = "/home/liudi/Dropbox/LocalTests/OutputMetricGraphRec/tables.csv";
    std::string folder = "/home/liudi/Dropbox/LocalTests/TwoAlgorithmTest/Outputs/";
//! We initilize a star:
    SingleStar star8(M_PI/6,8,4000,5,100,1,"Star8");
    SingleStar star4(M_PI/3,4,2000,5,100,1,"Star4");

//! We initilize filewriter:
    PrintToFile printer(folder);

// int proportion, double densityE, double parC
    MetricRec rec(2);
    rec.setTimePrecision(3);
    rec.addPostRunner(&printer);


    FastMSTLauncher mstComputator;
    mstComputator.setTimePrecision(2);
    mstComputator.addPostRunner(&printer);


    controller control(qq);

    control.addAlgorithm(&rec);
    control.addAlgorithm(&mstComputator);

    control.addCloud(&star4);
   ClassicDistanceMeasure distanceMeasure(2);
    CorrectEndTypeMeasure endTypeMeasure(2);
    NumberOfVertexMeasure vertexMeasure(2);
    CorrectTypeMeasure TypeMeasure(2);
    DualTreeDistanceMeasure newDistanceMeasure(2);
    MeanSquareDistanceMeasure ndm(2);

    control.addMeasurer(distanceMeasure);
    control.addMeasurer(endTypeMeasure);
    control.addMeasurer(vertexMeasure);
    control.addMeasurer(TypeMeasure);
    // control.addMeasurer(ndm);
    // control.addMeasurer(newDistanceMeasure);
    control.BeginTestRun();




}

void NewCoverTreeTree()
{
    using namespace mlpack;
    using namespace mlpack::neighbor; // NeighborSearch and NearestNeighborSort
    using namespace mlpack::metric;
    // Load the data from data.csv (hard-coded).  Use CLI for simple command-line
    // parameter handling.
    arma::mat referencedata;
    std::vector<Segment> segments;
    for (int i = 0; i < 10000; i++)
    {
        segments.push_back(Segment(Point(0,0,i),Point(0,0,i+1)));
    }

    GeneralConvertor::SegmentsToMat(segments,referencedata);
    arma::mat querydata;
    std::vector<Segment> segmentPoints;
    for (int j = 0; j < 10001; j++)
    {
        segmentPoints.push_back(Segment(Point(0,0,j),Point(0,0,j)));
    }

    GeneralConvertor::SegmentsToMat(segmentPoints,querydata);
//   data::Load("data.csv", data, true);
    // Use templates to specify that we want a NeighborSearch object which uses
    // the Manhattan distance.
    //  NeighborSearch<NearestNeighborSort, ManhattanDistance> nn(data);
    mlpack::neighbor::NeighborSearch<NearestNeighborSort,SegmentDistance,arma::mat,mlpack::tree::StandardCoverTree > nn(referencedata);

    // Create the object we will store the nearest neighbors in.
    arma::Mat<size_t> neighbors;
    arma::mat distances; // We need to store the distance too.
    // Compute the neighbors.
    nn.Search(querydata,1, neighbors, distances);
    // Write each neighbor and distance using Log.
    for (size_t i = 0; i < neighbors.n_elem; ++i)
    {
        if (distances[i] != 0)
        {
            std::cout << "Nearest neighbor of point " << i << " is segment"
                      << neighbors[i] << " and the distance is " << distances[i] << ".\n";
        }
    }
    //  std::cout << "Distancedebug: " << sqrt(CGAL::squared_distance(Segment[2],Point(0,10,10))) << std::endl;
//   double d = sqrt(CGAL::squared_distance(segments[2],Point(0,10,10))) ;
//   std::cout << d << std::endl;
//  std::cout << "Extra point" << std::endl;
}


void DebugMetricRec()
{
    std::string pathin = "/home/liudi/Dropbox/Skeletonization Project/OldFolders/unileverdata/XYZ_Files/Cluster_Frame00000002.xyz";
  // std::string pathin = "/home/liudi/Downloads/DataSkeletonization-master/out.xyz";
    std::list<Point> points;
    GeneralConvertor::XYZtoPoint(pathin,points);
    MetricRec rec(2);
    MyGraphType G;
    rec.Run(points,G);
    std::string pathout = "/home/liudi/Dropbox/LocalTests/OutputMetricGraphRec/out1.vtk";
    GeneralConvertor::GraphToVtk(pathout,G);





}

int main()
{
//Tests8();
//Test9();
//MSTTest();
//AMSTTest();
//AMSTTest();

//! This is required, to get proper random number sequence
//   srand( time( NULL ) );
    //  ControllerTest();

//NewCoverTreeTree();
///! We use this number sequence to debug the code:
//srand(128);
//std::cout <<rand()::numeric_limit<unit>::min(); << std::endl;ComputeDeluanayTriangulation(MyGraphType & G, std::list<Point> & Vector)
//! Here we test:ControllerTest()
//runGradientDescendTester();
//RunSeriousTests();
//RunCyclicTests();
//TestNewAMSTTreeType();
//  MassiveConvertion();
// TestNewAMSTTreeType();
//    TestAbstractThings();
//MapperWorkingTest();
//PrecisionTest();
//   MlPackTimerTest();
//ControllerTest();
//ControllerTest();
    // std::cout << "Compilation succeful" << std::endl;
    // std::cout << "one more thing" << std::endl;
    //  std::cout << "Bug fixed, actually" << std::endl;
    //  FileSystemTest();
// ControllerTestRealDataSimple();
//  ControllerTest();
    // std::cout << "Succeful compilation xD" << std::endl;
//ControllerTestRealDataSimple();

    // ControllerTest();
//   TryOutDeluanay();
    //  TryOutDeluanayCyclicTriangle();
MetriccontrollerTest();
 //   DebugMetricRec();
    std::cout << "Mlpackversion: " << mlpack::util::GetVersion() << std::endl;
    std::cout << "Compilation succeful" << std::endl;
    //  std::cout << "Compilation succeful" << std::endl;
//  PrintingStarCloudsTest();
//  PrintingStarCloudsTest();
    //  DistanceBetweenLineSegments();
//  ControllerTestRealDataSimple();
}

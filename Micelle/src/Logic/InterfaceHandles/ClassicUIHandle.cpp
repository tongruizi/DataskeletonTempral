#include "ClassicUIHandle.h"
#include "Filereader.h"
#include "Computation.h"
#include "FlexibleComplex.h"
#include "BranchDetection.h"
#include "StraighteningMethods.h"
#include "BranchSimplification.h"
#include "GraphGeneration.h"
#include "CloudGenerator.h"

ClassicUIHandle::ClassicUIHandle(double a, double b, double c)
{
    ClassicUIHandle::branchDetect = a;
    ClassicUIHandle::Straightening = b;
    ClassicUIHandle::branchCollapse = c;
}

void ClassicUIHandle::Read(std::list<Point> & cloud, std::string xyzfile)
{
    Filereader p;
    p.setPath(xyzfile);
    p.XYZRead(xyzfile,cloud);
}

//void ClassicUIHandle::MSTAllocationRecognition(MyGraphType & G, std::list<Point> & cloud)
//{
//tree = Computation::computeMST(cloud);
//BranchDetection::SimplifyIt(tree,optiout,param1,"");
//
//
//}

void ClassicUIHandle::CorrectCalculationToGraphIntermediate(MyGraphType & FG, std::list<std::list<Point>> & optipath,
std::list<std::list<Point>> & branchsimplified, std::list<Point> & cloud, double & mstlength, std::string settings)
{
MyGraphType tree = Computation::computeMST(cloud);
MyGraphType optiout;
mstlength = Computation::AverageEdgelength(tree);
double param1 = ClassicUIHandle::branchDetect*mstlength;
BranchDetection::SimplifyIt(tree,optiout,param1,"",settings);
double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, this->Straightening );
std::cout << "The straightening parameter: " << ddd << std::endl;
double valuev = ddd*ClassicUIHandle::branchCollapse;
BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
BranchSimplification::PathToGraphProper(FG, branchsimplified);
}

void ClassicUIHandle::CorrectCalculationToGraphMSTsave(MyGraphType & FG, std::list<Point> & cloud, double & mstlength, MyGraphType & tree)
{
tree = Computation::computeMST(cloud);
MyGraphType optiout;
mstlength = Computation::AverageEdgelength(tree);
double param1 = ClassicUIHandle::branchDetect*mstlength;
BranchDetection::SimplifyIt(tree,optiout,param1,"","pure");
std::list<std::list<Point>> optipath;
std::list<std::list<Point>> branchsimplified;
double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, ClassicUIHandle::Straightening);
double valuev = ddd*ClassicUIHandle::branchCollapse;
BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
BranchSimplification::PathToGraphProper(FG, branchsimplified);
}

void ClassicUIHandle::CorrectCalculationToGraph(MyGraphType & FG, std::list<Point> & cloud, double & mstlength, std::string settings)
{
MyGraphType tree = Computation::computeMST(cloud);
MyGraphType optiout;
mstlength = Computation::AverageEdgelength(tree);
double param1 = ClassicUIHandle::branchDetect*mstlength;
BranchDetection::SimplifyIt(tree,optiout,param1,"",settings);
std::list<std::list<Point>> optipath;
std::list<std::list<Point>> branchsimplified;
double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, this->Straightening );
std::cout << "The straightening parameter: " << ddd << std::endl;
double valuev = ddd*ClassicUIHandle::branchCollapse;
BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
BranchSimplification::PathToGraphProper(FG, branchsimplified);
}

void ClassicUIHandle::CorrectCalculation(std::string & in, std::string & out, double & mstlength)
{
std::list<Point> cloud;
ClassicUIHandle::Read(cloud, in);
std::cout << "Cloud size: " << cloud.size() << std::endl;
MyGraphType tree = Computation::computeMST(cloud);
MyGraphType optiout;
mstlength = Computation::AverageEdgelength(tree);
double param1 = ClassicUIHandle::branchDetect*mstlength;
BranchDetection::SimplifyIt(tree,optiout,param1,"","pure");
std::list<std::list<Point>> optipath;
std::list<std::list<Point>> branchsimplified;

double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, ClassicUIHandle::Straightening);
double valuev = ddd*ClassicUIHandle::branchCollapse;
BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
Write::pathPrintToVtkPointlist(branchsimplified,out);

}

void ClassicUIHandle::CoreCalculation(std::list<Point> & cloud, std::string f1, std::string f2,
                                      std::string f3, std::string f4, std::string f5, std::string f6)
{
    std::cout << "Paramters:" << std::endl;
    std::cout << "Branch detect: " << ClassicUIHandle::branchDetect << std::endl;
    std::cout << "Straightening parameter: " << ClassicUIHandle::Straightening << std::endl;
    std::cout << "Branch parameter: " << ClassicUIHandle::branchCollapse << std::endl;
    Alpha_shape_3 as(cloud.begin(),cloud.end());
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
    MyGraphType tree = Computation::computeMST(cloud);
    Write::GraphToVtk(f1,tree);
    Write::AlphaVTK(f2, as);
//MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as, f3);
//Write::GraphToVtk(f3,ABC);
    MyGraphType optiout;
    BranchDetection::SimplifyIt(tree,optiout,ClassicUIHandle::branchDetect,"","pure");
    Write::GraphToVtk(f4,optiout);
    std::list<std::list<Point>> optipath;
    double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, ClassicUIHandle::Straightening);
    Write::pathPrintToVtkPointlist(optipath,f5);
    std::list<std::list<Point>> branchsimplified;
    double valuev = ddd*ClassicUIHandle::branchCollapse;
    BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
    Write::pathPrintToVtkPointlist(branchsimplified,f6);

// Calculating stuff here
}

void ClassicUIHandle::CalculateAll(std::string in, std::string f1, std::string f2,
                                   std::string f3, std::string f4, std::string f5, std::string f6)
{
    std::list<Point> cloud;
    ClassicUIHandle::Read(cloud, in);
    ClassicUIHandle::CoreCalculation(cloud,f1,f2,f3,f4,f5,f6);
//StraighteningMethods::Optimize(optipath, StraighteningMethods::GraphToPaths(pathsvariable))
}

void ClassicUIHandle::runSingleTest(int n, int k, double minangle, double scale, double epsilon,std::vector<std::string> & files)
{
    std::list<Point> out;
    MyGraphType G;
    GraphGeneration::RandomGraph1(k,minangle,scale,G);
    CloudGenerator::generatePoints(n, G, epsilon, out);
    ClassicUIHandle::CoreCalculation(out, files[0], files[1], files[2], files[3], files[4], files[5]);
}

void ClassicUIHandle::multipleTests(int n, int k, double minangle, double scale, int testN, double epsilon,std::string folder)
{
    for (int i = 0; i < testN; i++)
    {
        std::vector<std::string> files;
        std::string f1 = folder + "mst" + std::to_string(i) + ".vtk";
        std::string f2 = folder + "alpha" + std::to_string(i) + ".vtk";
        std::string f3 = folder + "optimst" + std::to_string(i) + ".vtk";
        std::string f4 = folder + "optipath" + std::to_string(i) + ".vtk";
        std::string f5 = folder + "simplified" + std::to_string(i) + ".vtk";
        std::string f6 = folder + "branchsimplified" + std::to_string(i) + ".vtk";
        files.push_back(f1);
        files.push_back(f2);
        files.push_back(f3);
        files.push_back(f4);
        files.push_back(f5);
        files.push_back(f6);
        runSingleTest(n,k,minangle,scale,epsilon,files);
    }
}


bool ClassicUIHandle::CorrectFormSecond(MyGraphType & G, int n, int m)
{
    std::set<int> setInt;
    int n2 = n + 1;
    int m2 = m + 1;
    setInt.insert(n2);
    setInt.insert(m2);
    int onev;
    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        int tmpD = boost::degree(*it, G);
        if (tmpD > 2)
        {
            std::cout << "Currently: " << tmpD << std::endl;
            auto bt = setInt.find(tmpD);
            if (bt == setInt.end())
            {
            std::cout << tmpD << std::endl;
            std::cout << "EXIT HERE" << std::endl;
                return false;
            }
            else
            {
            std::cout << "ERASED SOMETHING" << std::endl;
            setInt.erase(tmpD);
            }
        }
        else if (tmpD == 1)
        {
            onev++;

        }

    }
    std::cout << "SETINT :" << setInt.size() << " | " << "ONEV:" << onev << std::endl;
    return (setInt.size() == 0) && (onev == n + m);


}

bool ClassicUIHandle::CorrectFormFirst(MyGraphType & G, int n)
{
    std::cout << "Parameter N:" << n << std::endl;
    auto vpair = boost::vertices(G);
    int degree1 = 0;
    int starn = 0;
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        int tmpD = boost::degree(*it, G);
        if (tmpD == 1)
        {
            degree1++;
        }
        else if (tmpD == n)
        {
            starn++;
        }
        else if ((tmpD != 2) && (tmpD != n) && (tmpD != 1))
        {
            std::cout << "Unsucceful, tmpD =" << tmpD << std::endl;
            return false;
        }
    }
    std::cout << "Stats: N-star: " << starn << " | 1-vertex: " << degree1 << std::endl;
    return (starn == 1) && (degree1 == n);

}


void ClassicUIHandle::GraphRecognitionTestFirst(int n, int k, double minangle, double scale, int testN, double epsilon, std::string folder)
{
    int correct = 0;
    for (int i = 0; i < testN; i++)
    {
        std::cout << "Round " << i << " beginning" << std::endl;
        std::list<Point> cloud;
        std::string f1 = folder + "cloud" + std::to_string(i) + ".vtk";
        MyGraphType G;
        GraphGeneration::RandomGraph1(k,minangle,scale,G);
        CloudGenerator::generatePoints(n, G, epsilon, cloud);
        MyGraphType tree = Computation::computeMST(cloud);
        //Print the tree:
        Write::GraphToVtk(f1,tree);

        Alpha_shape_3 as(cloud.begin(),cloud.end());
        Alpha_iterator opt = as.find_optimal_alpha(1);
        as.set_alpha(*opt);
        MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as, "");
        MyGraphType optiout;
        std::cout << "The parameter: " << ClassicUIHandle::branchDetect << std::endl;
        BranchDetection::SimplifyIt(ABC,optiout,ClassicUIHandle::branchDetect,"","pure");
        std::list<std::list<Point>> optipath;
        std::string f3 = folder + "longestpath" + std::to_string(i) + ".vtk";
        std::string f4 = folder + "optipath" + std::to_string(i) + ".vtk";
        Write::GraphToVtk(f3,optiout);
        double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, ClassicUIHandle::Straightening);
        double valuev = ddd*ClassicUIHandle::branchCollapse;
        std::list<std::list<Point>> branchsimplified;
        BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
        std::string f2 = folder + "final" + std::to_string(i) + ".vtk";
        // Print the output:
        Write::pathPrintToVtkPointlist(branchsimplified,f2);

        MyGraphType FG;
        BranchSimplification::PathToGraph(FG, branchsimplified);
        if (ClassicUIHandle::CorrectFormFirst(FG,k))
        {
            correct++;

        }
    }
    std::cout << "Number of correct recognitions: " << correct << std::endl;

}

void ClassicUIHandle::GraphRecognitionTestSecond(int n, int k, int k2, double minangle, double scale, int testN, double epsilon, std::string folder)
{
    int correct = 0;
    for (int i = 0; i < testN; i++)
    {
        std::cout << "Round " << i << " beginning" << std::endl;
        std::list<Point> cloud;
        std::string f1 = folder + "cloud" + std::to_string(i) + ".vtk";
        MyGraphType G;
        GraphGeneration::RandomGraph2(k,k2,minangle,scale,G);
        CloudGenerator::generatePoints(n, G, epsilon, cloud);
        MyGraphType tree = Computation::computeMST(cloud);
        //Print the tree:
        Write::GraphToVtk(f1,tree);

        Alpha_shape_3 as(cloud.begin(),cloud.end());
        Alpha_iterator opt = as.find_optimal_alpha(1);
        as.set_alpha(*opt);
        MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as, "");
        MyGraphType optiout;
        std::cout << "The parameter: " << ClassicUIHandle::branchDetect << std::endl;
        BranchDetection::SimplifyIt(ABC,optiout,ClassicUIHandle::branchDetect,"","pure");
        std::list<std::list<Point>> optipath;
        std::string f3 = folder + "longestpath" + std::to_string(i) + ".vtk";
        std::string f4 = folder + "optipath" + std::to_string(i) + ".vtk";
        Write::GraphToVtk(f3,optiout);
        double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, ClassicUIHandle::Straightening);
        double valuev = ddd*ClassicUIHandle::branchCollapse;
        Write::pathPrintToVtkPointlist(optipath,f4);
        std::list<std::list<Point>> branchsimplified;
        BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
        std::string f2 = folder + "final" + std::to_string(i) + ".vtk";
        // Print the output:
        Write::pathPrintToVtkPointlist(branchsimplified,f2);

        MyGraphType FG;
        BranchSimplification::PathToGraph(FG, branchsimplified);
        if (CorrectFormSecond(FG, k, k2))
        {
            std::cout << "Correct" << std::endl;
            correct++;

        }
        else
        {
        std::cout << "INCORRECT" << std::endl;
        }
    }
    std::cout << "Number of correct recognitions: " << correct << std::endl;

}




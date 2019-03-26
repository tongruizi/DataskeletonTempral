#ifndef UIHANDLER_H
#define UIHANDLER_H

#include "Definitions.h"
#include "Mapper_Parameters.h"
#include "AlphaReeb_Parameters.h"
#include "Graph.h"

class UIHandler
{
    std::string folder;
    public:
        UIHandler();
        void setFolder(std::string folder);
        void multipleTests(int n, int k, double minangle, double scale, int testN, double epsilon);
        void runTest(int n, int k, double minangle, double scale, std::string filename, double epsilon, std::string f2, std::string f3, std::string f4);
        static void AlphaMST(std::list<Point> cloud, std::string f1, std::string f2, std::string f3,std::string f4);
        static void AlphaSimplification(std::list<Point> cloud, std::string f1, std::string f2);
        void runTestPremium(int n, int k, double minangle, double scale, std::string filename, double epsilon, std::string f2, std::string fi);
        void multipleTestsPremium(int n, int k, double minangle, double scale, int testN, double epsilon);
        static void AlphaDistanceComputer(std::list<Point> cloud, std::string f1, std::string f2, std::string f3, std::string f4);
        static void HashComplexTest(std::list<Point> cloud);
        static void AlphaReebTest(std::list<Point> cloud, std::string outwrite, std::string inter);
        static void MapperTest(std::string inwrite, std::string outwrite, Mapper_Parameters const & param);
        static void ComputeAlphaReeb(std::string inwrite, std::string outwrite, AlphaReeb_Parameters & param, double epsilon);
        static void ThreeAlgorithmSuperTest(std::string inputfolder, std::string outputfolder, std::vector<std::string> inputs, std::vector<std::string> outputs);
        static void ComputeAlphaReebInterior(MyGraphType & out, AlphaReeb_Parameters & param, std::vector<Point> & cloud, double epsilon);
        static void RunTestsForTheeAlgorithm(int n, int k, double minangle, double scale, int testN, double epsilon);
        static void DebugClassic(int n, int k, double minangle, double scale, int testN, double epsilon);
        static void BetterTesting(int n, int k, double minangle, double scale, int testN, double epsilon);
        static void FindBug(int n, int k, double minangle, double scale, int testN, double epsilon);
        static void RunTheFinalTests(int k, double minangle, double scale, int testN, double epsilon, double alpha, double mappercluster);
        static void ComputeThingsRequiredForPaper();

       // static void ComputeMapper(std::string inwrite, std::string outwrite);

    protected:

    private:
};

#endif // UIHANDLER_H

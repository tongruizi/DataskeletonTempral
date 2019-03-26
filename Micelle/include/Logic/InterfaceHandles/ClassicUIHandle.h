#ifndef CLASSICUIHANDLE_H
#define CLASSICUIHANDLE_H
#include "Definitions.h"
#include "Write.h"

class ClassicUIHandle
{
    double branchDetect;
    double Straightening;
    double branchCollapse;

    public:
        ClassicUIHandle(double, double, double);
        void Read(std::list<Point> & cloud, std::string path);
        void CalculateAll(std::string in, std::string f1, std::string f2,
        std::string f3, std::string f4, std::string f5, std::string f6);
        void multipleTests(int n, int k, double minangle, double scale, int testN, double epsilon,std::string folder);
        void runSingleTest(int n, int k, double minangle, double scale, double epsilon,std::vector<std::string> & files);
        void CoreCalculation(std::list<Point> & cloud, std::string f1, std::string f2,
          std::string f3, std::string f4, std::string f5, std::string f6);
        static bool CorrectFormFirst(MyGraphType & G, int n);
        void GraphRecognitionTestFirst(int n, int k, double minangle, double scale, int testN, double epsilon, std::string folder);
        void GraphRecognitionTestSecond(int n, int k, int k2, double minangle, double scale, int testN, double epsilon, std::string folder);
        static bool CorrectFormSecond(MyGraphType & G, int n, int m);
        void CorrectCalculation(std::string & in, std::string & out, double & mstl);
        void CorrectCalculationToGraph(MyGraphType & FG, std::list<Point> & cloud, double & mstlength,std::string settings);
        void CorrectCalculationToGraphMSTsave(MyGraphType & FG, std::list<Point> & cloud, double & mstlength, MyGraphType & mst);
        void CorrectCalculationToGraphIntermediate(MyGraphType & FG, std::list<std::list<Point>> & optipath, std::list<std::list<Point>> & branchsimplified
               , std::list<Point> & cloud, double & mstlength, std::string settings);


      //  smethod.CorrectCalculationToGraph(outGraph, parameters[j], cloud, mstlength,"sd");

    protected:

    private:
};

#endif // CLASSICUIHANDLE_H

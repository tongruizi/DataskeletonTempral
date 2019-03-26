#ifndef ABSTRACTALGORITHM_H
#define ABSTRACTALGORITHM_H

#include "Definitions.h"
#include "Graph.h"
#include "ExplicitMeasurer.h"
#include <chrono>
#include "PostRunInterface.h"
#include "ManualMeasure.h"

class AbstractAlgorithm
{
public:
    //! Constructor:
    AbstractAlgorithm(std::string name):name(name),TimeMeasure("Time",1),RunNumber(0) {}
    //! Destructor:
    ~AbstractAlgorithm()
    {
        for (int i = 0; i < measurers.size(); i++)
        {
            delete measurers[i];
        }
    }

    virtual void Run(std::list<Point> & cloudlist, MyGraphType & out) = 0;

    void ExplicitRun(std::list<Point> & cloudlist, MyGraphType & out, AbstractCloudType* gen)
    {
        //! Something will be add here soon.
        auto start = std::chrono::system_clock::now();
        //! Algorithm itself:
        Run(cloudlist, out);
        //! Time counting methods:
        auto endd = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = endd-start;
        TimeMeasure.CustomAppend(diff.count());
        //! PostRuners:
        for (int j = 0; j < postRunners.size(); j++)
        {
            //run(MyGraphType & G, generatable* gen, int RunNumber, std::string AlgorithmName);
            postRunners[j] -> run(out,cloudlist,gen,RunNumber,name);
        }
        //! Update Run number:
        RunNumber++;


    }

    std::string returnName()
    {
        return name;
    }
    int returnMeasurerSize()
    {
        return measurers.size();
    }
    void addMeasurer(ExplicitMeasurer* m)
    {
        measurers.push_back(m);
    }
    ExplicitMeasurer* returnMeasurer(int k)
    {
        return measurers[k];
    }
    void CycleOverMeasurers(MyGraphType & G, AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)
    {
        for (int k = 0; k < measurers.size(); k++)
        {
          //  std::cout << "Measurment "<< k << " initilized succefully...." << std::endl;
            (*measurers[k]).run(G,cloud,generatedCloud,cloudit);
          //  std::cout << "Measurment "<< k << " Completed succefully!!!" << std::endl;
        }
    }
    void ResetAllMeasurers()
    {
        for (int k = 0; k < measurers.size(); k++)
        {
            measurers[k] ->resetStatistic();
        }
        TimeMeasure.resetStatistic();
    }

    void setRunNumber(int k)
    {
        RunNumber = k;
    }

    void WriteResults(std::vector<std::string> & rvector)
    {
        rvector.resize(measurers.size());
        for (int k = 0; k < measurers.size(); k++)
        {
            rvector[k] = measurers[k] ->returnStatisticString() ;
        }
    }
    std::string returnTimeMeasureString()
    {
        return TimeMeasure.returnStatisticString();
    }

    void addPostRunner(PostRunInterface* postRunner)
    {
        postRunners.push_back(postRunner);
    }

    void setTimePrecision(int p)
    {
        TimeMeasure.setPrecision(p);
    }



protected:
    std::vector<ExplicitMeasurer*> measurers;
    ManualMeasure<double> TimeMeasure;
    //! We will create also method for post processing, such as printing
    std::vector<PostRunInterface*> postRunners;
private:
    std::string name;
    int RunNumber;
};

#endif // ABSTRACTALGORITHM_H

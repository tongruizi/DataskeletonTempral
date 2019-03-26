#include "controller.h"
#include "generatable.h"
#include <thread>
#include "Definitions.h"
#include "GeneralConvertor.h"
//#define THREADS 7;

controller::controller(std::string filename):
    filename(filename)
{
}

void controller::InitilizeResultsCube()
{
    results = std::vector<std::vector<std::vector<std::string>>>(ACT.size(), std::vector<std::vector<std::string>>(algorithms.size(), std::vector<std::string>(measurerNames.size(),"")));
    runTimeStorage = std::vector<std::vector<std::string>>(ACT.size(),std::vector<std::string>(algorithms.size(),""));
}

void controller::WriteDownToLatexTable()
{

    std::vector<std::string> cloudNames(ACT.size());
    std::vector<std::string> algorithmNames(algorithms.size());

//! Convert them here:

    for (int i = 0; i < ACT.size(); i++)
    {
        cloudNames[i] = ACT[i] -> returnName();
    }
    for (int i = 0; i < algorithms.size(); i++)
    {
        algorithmNames[i] = algorithms[i] ->returnName();
    }
//! Just to make everything work, this is done in non-clean way.

    GeneralConvertor::DataToCSV(results,runTimeStorage,cloudNames,algorithmNames,measurerNames,filename);

}


void controller::InsideCloudListLoop(int i, int clit)
{
    std::list<Point> points;
    //! Measure the time!
    auto start = std::chrono::system_clock::now();
    //! Algorithm itself:
    ACT[i]->GenerateCloud(points,clit);
    //! Time counting methods:
    auto endd = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = endd-start;
    std::cout << "Time on generating the cloud: " << diff.count() << std::endl;

    for (int j = 0; j < algorithms.size(); j++)
    {
        //! Print where are we right now:
        std::cout << "Dealing with: ";
        std::cout << " CloudType: ";
        std::cout << std::to_string(i);
        std::cout << " Algorithm: " << std::to_string(j) << " IterationNumber: " << std::to_string(clit) << std::endl;
        //! Step 2: Run Algorithm
        MyGraphType out;
        algorithms[j]->ExplicitRun(points,out,ACT[i]);

        //! Step 3: Do measurments:

        algorithms[j]->CycleOverMeasurers(out,ACT[i], &points,clit);
    }

}

void controller::LoopOverCloudList(int i)
{
//! Section below will be rewritten for multi threading application. Idea would be that we would allocate

    for (int clit = 0; clit < (*ACT[i]).returnNumberOfRuns(); clit++ )
    {

        //i presume this method can be used here
        //! The inside function is written as separate instance, so we can launch threads involving this
        this->InsideCloudListLoop(i,clit);



    }
    //! If this will be threading application, then we will wait here for all proccesses to halt and then write down the results.
    //! Step 4 Write down results. Reset each measurer
    for (int j = 0; j < algorithms.size(); j++)
    {
        std::vector<std::string> rvector;
        algorithms[j] -> WriteResults(rvector);
        results[i][j] = rvector;
        runTimeStorage[i][j] = algorithms[j] -> returnTimeMeasureString();


        //! Here we reset stuff
        algorithms[j] -> ResetAllMeasurers();
        algorithms[j] -> setRunNumber(0);
    }

}

void controller::BeginTestRun()
{
    this->InitilizeResultsCube();
    for (int i = 0; i < ACT.size(); i++)
    {
        //this->InitilizeAllMeasurments();

        this->LoopOverCloudList(i);
    }
    //! We will write everything into files
    this->WriteDownToLatexTable();
    std::cout << "Write succeful" << std::endl;

}

void controller::addAlgorithm(AbstractAlgorithm* k)
{
    algorithms.push_back(k);
}
void controller::addCloud(AbstractCloudType* k)
{
    ACT.push_back(k);
}
void controller::addMeasurer(ExplicitMeasurer & q)
{
    for (int i = 0; i < algorithms.size(); i++)
    {
        algorithms[i] ->addMeasurer(q.Clone());
    }
    measurerNames.push_back(q.returnName());

}

void controller::FlushClouds(std::vector<CloudTypePrinterAlgorithm*> & p)
{
    for (int i = 0; i < this->ACT.size(); i++)
    {
        p[i] ->CompleteOperation(ACT[i]);
    }



}


//
//void controller::LoopOverMeasurers(MyGraphType & G, int i, int j, std::list<Point> & cloud)
//{
//
//    //! set size of Information vector
//    for (int k = 0; k < measurers.size(); k++)
//    {
//
//      (*measurers[k]).run(j, G ,ACT[i],&cloud);
//
//    }
//
//}

//void controller::InitilizeAllMeasurments()
//{
//
//    for (int k = 0; k < measurers.size(); k++)
//    {
//    measurers[k] ->resetStatistic(algorithms.size());
//    }
//
////    for (auto MR = measurers.begin(); MR != measurers.end(); MR++)
////    {
////        (*MR)->resetStatistic(algorithms.size());
////        // (*MR)->setCloud(cloud);
////        //   std::string returnStatisticData();
////    }
//
//}






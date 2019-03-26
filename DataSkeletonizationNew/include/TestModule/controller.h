#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "AbstractAlgorithm.h"
#include "generatable.h"
#include "AbstractMeasurer.h"
#include "ExplicitMeasurer.h"
#include "AbstractCloudType.h"
#include "CloudTypePrinterAlgorithm.h"

//void controller::InsideCloudListLoop(int i, int clit)

//void controller::LoopOverCloudList(int i)

class controller
{
public:
    controller(std::string filename);
    void InitilizeResultsCube();
    void WriteDownToLatexTable();
  //  void LoopOverMeasurers(MyGraphType & G, int i, int j, std::list<Point> & cloud);
  //  void InitilizeAllMeasurments();
    void InsideCloudListLoop(int i, int clit);
    void LoopOverCloudList(int i);
    void BeginTestRun();

    void addAlgorithm(AbstractAlgorithm* k);
    void addCloud(AbstractCloudType* k);
    void addMeasurer(ExplicitMeasurer & q);

    void FlushClouds(std::vector<CloudTypePrinterAlgorithm*> & p);
    //void
protected:

private:
    int msize;
    std::string filename;
    //! Vector for cloudtypes:
    std::vector<AbstractCloudType*> ACT;
    //! Vector for algorithmtypes:
    std::vector<AbstractAlgorithm*> algorithms;

    //! Table of results in std::string format
    std::vector<std::vector<std::vector<std::string>>> results;
    std::vector<std::vector<std::string>> runTimeStorage;
    //! Coordinate (i, j, k) represents Cloudtype i, Algorithm j and results from measure k.

    //! Vector consisting of name measurers.
    std::vector<std::string>  measurerNames;

   // StatisticInfo<int> correctRecognition;
   // StatisticInfo<double> MinMaxError;
   // std::vector<Measurer> measures;
};

#endif // CONTROLLER_H

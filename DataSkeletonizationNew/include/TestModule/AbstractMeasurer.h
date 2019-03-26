#ifndef ABSTRACTMEASURER_H
#define ABSTRACTMEASURER_H

#include "StatisticInfo.h"
#include "ExplicitMeasurer.h"
#include "AbstractCloudType.h"

template <class T>
class AbstractMeasurer : public ExplicitMeasurer
{
public:
    AbstractMeasurer(std::string name,int precision):ExplicitMeasurer(name,precision)
    {

    }

    virtual T CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit) =0;
    virtual std::string returnStatisticString() = 0;

    void run(MyGraphType & G, AbstractCloudType* cloud,  std::list<Point>* generatedCloud, int cloudit) override
    {
        T InMeasure=CompleteMeasurments(G,cloud,generatedCloud,cloudit);
        mystatistic.add(InMeasure);
    }
    //! No longer needed
//    T getStatistic_instance()
//    {
//        return mystatistic;
//    }
    std::string returnStatisticData()
    {
        return mystatistic.returnInfo();
    }

    void resetStatistic() override
    {
        mystatistic.reset();
    }
    StatisticInfo<T>* provideAccessToMyStatistic()
    {
    return &mystatistic;
    }

    virtual AbstractMeasurer<T>* Clone() = 0;

protected:
    StatisticInfo<T> mystatistic;
private:


};




#endif // ABSTRACTMEASURER_H

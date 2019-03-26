#ifndef EXPLICITMEASURER_H
#define EXPLICITMEASURER_H

#include "Graph.h"
#include "generatable.h"
#include "AbstractCloudType.h"

//! this is the true interface class
class ExplicitMeasurer
{
public:
    ExplicitMeasurer(std::string name,int precision):
        name(name),precision(precision) {}
    virtual std::string returnStatisticString() = 0;
    virtual void run(MyGraphType & G, AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)=0;
    virtual void resetStatistic() = 0;
    virtual ExplicitMeasurer* Clone() = 0;
//    void setCloud(generatable* cloudt)
//    {
//    cloudt = cloud;
//    }
//    void setGeneratedCloud(std::list<Point>*  cloudt)
//    {
//    generatedCloud = cloudt;
//    }
    std::string returnName()
    {
        return name;
    }
    void setPrecision(int q)
    {
        precision = q;

    }

    // virtual double SimplifiedOutput() = 0;
protected:
    std::string name;
    int precision;
private:

};

#endif // EXPLICITMEASURER_H

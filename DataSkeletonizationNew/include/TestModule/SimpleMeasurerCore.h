#ifndef SIMPLEMEASURERCORE_H
#define SIMPLEMEASURERCORE_H

#include "SimpleStatistic.h"
#include "ExplicitMeasurer.h"

class SimpleMeasurerCore : public ExplicitMeasurer
{
public:
    SimpleMeasurerCore(std::string name, int precision):
        ExplicitMeasurer(name,precision) {}


    std::string returnStatisticString() override
    {

        double avg = statistic.returnAvg();
        avg = 100 * avg;
        std::stringstream stream;
        stream << std::fixed << std::setprecision(precision) << avg;
        std::string s = stream.str() + "%";
        return s;
    }

    virtual bool CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit) =0;

    void run(MyGraphType & G, AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit) override
    {
        bool ww = this->CompleteMeasurments(G,cloud,generatedCloud,cloudit);
        statistic.append(ww);
    }


    void resetStatistic() override
    {
        statistic.reset();
    }
    virtual SimpleMeasurerCore* Clone() = 0;

protected:
    SimpleStatistic statistic;

private:

};

#endif // SIMPLEMEASURERCORE_H

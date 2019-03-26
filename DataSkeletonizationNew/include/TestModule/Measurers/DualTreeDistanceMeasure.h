#ifndef DUALTREEDISTANCEMEASURE_H
#define DUALTREEDISTANCEMEASURE_H

#include "AbstractMeasurer.h"


class DualTreeDistanceMeasure : public AbstractMeasurer<double>
{
public:
    DualTreeDistanceMeasure(int precision);
    double CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit);
    std::string returnStatisticString();
    DualTreeDistanceMeasure* Clone() override
    {
        return new DualTreeDistanceMeasure(*this);
    }

protected:

private:
};

#endif // DUALTREEDISTANCEMEASURE_H

#ifndef MEANSQUAREDISTANCEMEASURE_H
#define MEANSQUAREDISTANCEMEASURE_H

#include "AbstractMeasurer.h"

class MeanSquareDistanceMeasure : public AbstractMeasurer<double>
{
public:
    MeanSquareDistanceMeasure(int precision);
    double CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit);
    std::string returnStatisticString();
    MeanSquareDistanceMeasure* Clone() override
    {
        return new MeanSquareDistanceMeasure(*this);
    }


protected:

private:
};

#endif // MEANSQUAREDISTANCEMEASURE_H

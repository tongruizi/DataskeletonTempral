#ifndef CLASSICDISTANCEMEASURE_H
#define CLASSICDISTANCEMEASURE_H

#include "AbstractMeasurer.h"
#include "AbstractCloudType.h"

class ClassicDistanceMeasure : public AbstractMeasurer<double>
{
public:
    ClassicDistanceMeasure(int precision);
    double CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit);
    std::string returnStatisticString();
    ClassicDistanceMeasure* Clone() override
    {
        return new ClassicDistanceMeasure(*this);
    }



protected:

private:
};

#endif // CLASSICDISTANCEMEASURE_H

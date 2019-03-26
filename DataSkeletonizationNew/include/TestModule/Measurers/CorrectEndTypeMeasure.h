#ifndef CORRECTENDTYPEMEASURE_H
#define CORRECTENDTYPEMEASURE_H

#include "Graph.h"
#include "ExplicitMeasurer.h"
#include "SimpleStatistic.h"
#include "SimpleMeasurerCore.h"
#include "AbstractCloudType.h"

class CorrectEndTypeMeasure : public SimpleMeasurerCore
{
    public:
        CorrectEndTypeMeasure(int precision);
        bool CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit);
        CorrectEndTypeMeasure* Clone() override
        {
        return new CorrectEndTypeMeasure(*this);
        }
    protected:

    private:
};

#endif // CORRECTENDTYPEMEASURE_H

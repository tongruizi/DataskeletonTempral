#ifndef MANUALMEASURE_H
#define MANUALMEASURE_H

#include "AbstractMeasurer.h"
#include "AbstractCloudType.h"

template <class T>
class ManualMeasure : public AbstractMeasurer<T>
{
public:
    ManualMeasure(std::string name, int precision):
        AbstractMeasurer<T>(name,precision) {}
    virtual T CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)
    {
        //! We do nothing
        return 0;
    }
    std::string returnStatisticString()
    {
        double avg = (this -> mystatistic).returnAvg();
        std::stringstream stream;
        stream << std::fixed << std::setprecision(this->precision) << avg;
        std::string s = stream.str();
        return s;
    }
    void CustomAppend(T measurment)
    {
        (this -> provideAccessToMyStatistic()) -> add(measurment);
    }
    ManualMeasure<T>* Clone() override
    {
        return new ManualMeasure<T>(*this);
    }


protected:

private:
};

#endif // MANUALMEASURE_H

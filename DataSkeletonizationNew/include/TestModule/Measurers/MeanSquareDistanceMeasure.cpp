#include "MeanSquareDistanceMeasure.h"
#include "Computation.h"
#include <iostream>
#include <iomanip>

MeanSquareDistanceMeasure::MeanSquareDistanceMeasure(int precision):
    AbstractMeasurer("Distance",precision)
{
    //ctor
}

double MeanSquareDistanceMeasure::CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)
{
    return Computation::MeanSquareErrorAABB(G,*generatedCloud);
}

std::string MeanSquareDistanceMeasure::returnStatisticString()
{
    double avg = mystatistic.returnAvg();
    //  std::cout << "Returned: " << avg << std::endl;
    //  std::cout << "Data contained: " << mystatistic.returnNumber() << std::endl;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << avg;
    std::string s = stream.str();
    return s;
}

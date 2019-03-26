#include "ClassicDistanceMeasure.h"
#include "Computation.h"
#include <iostream>
#include <iomanip>

ClassicDistanceMeasure::ClassicDistanceMeasure(int precision):
    AbstractMeasurer("Distance",precision)
{
    //ctor
}

double ClassicDistanceMeasure::CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)
{
    double Error = Computation::AABBError(G,*generatedCloud);
//    std::cout << "Computed error: " << Error << std::endl;
    return Error;
}

std::string ClassicDistanceMeasure::returnStatisticString()
{
    double avg = mystatistic.returnAvg();
  //  std::cout << "Returned: " << avg << std::endl;
  //  std::cout << "Data contained: " << mystatistic.returnNumber() << std::endl;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << avg;
    std::string s = stream.str();
    return s;
}


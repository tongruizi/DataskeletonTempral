#include "CorrectEndTypeMeasure.h"
#include "Graph.h"


CorrectEndTypeMeasure::CorrectEndTypeMeasure(int precision):SimpleMeasurerCore("Endpoints",precision)
{
    //ctor
}


bool CorrectEndTypeMeasure::CompleteMeasurments(MyGraphType & G,AbstractCloudType* cloud, std::list<Point>* generatedCloud, int cloudit)
{
    bool theVariable = cloud->DoesGraphHaveCorrectForm(G,cloudit);
   // std::cout << "The descision (does graph have the correct form? " << theVariable << std::endl;
    return theVariable;
}






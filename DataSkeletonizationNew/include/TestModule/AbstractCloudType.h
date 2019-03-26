#ifndef ABSTRACTCLOUDTYPE_H
#define ABSTRACTCLOUDTYPE_H

#include "Graph.h"
#include "Definitions.h"

class AbstractCloudType
{
public:
    AbstractCloudType(int number_of_runs, std::string name):
        number_of_runs(number_of_runs), name(name)
    {}
    virtual bool IsGraphCorrect(MyGraphType & G, int iterationNumber) = 0;
    virtual bool DoesGraphHaveCorrectForm(MyGraphType & G, int iterationNumber) =0;
    virtual void GenerateCloud(std::list<Point> & p, int iterationNumber) = 0;

    int returnNumberOfRuns()
    {
        return number_of_runs;
    }
    std::string returnName()
    {
        return name;
    }
    void setNumberOfRuns(int k)
    {
    number_of_runs = k;
    }

protected:
    int number_of_runs;
    std::string name;


private:
};

#endif // ABSTRACTCLOUDTYPE_H

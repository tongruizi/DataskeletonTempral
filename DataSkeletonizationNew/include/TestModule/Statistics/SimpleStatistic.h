#ifndef SIMPLESTATISTIC_H
#define SIMPLESTATISTIC_H


class SimpleStatistic
{
public:
    SimpleStatistic():
        caseNumber(0), success(0), failure(0)
    {}

    void append(bool result)
    {
        caseNumber++;
        if (result)
        {
        success++;
        }
        else
        {
        failure++;
        }
    }
    int returnSuccess()
    {
    return success;
    }
    int returnFailure()
    {
    return failure;
    }
    int returnNumberOfRuns()
    {
    return caseNumber;
    }
    void reset()
    {
    caseNumber = 0;
    success = 0;
    failure = 0;
    }
    double returnAvg()
    {
    return (double) success/ (double) caseNumber;
    }



protected:

private:
    int caseNumber;
    int success;
    int failure;
};

#endif // SIMPLESTATISTIC_H

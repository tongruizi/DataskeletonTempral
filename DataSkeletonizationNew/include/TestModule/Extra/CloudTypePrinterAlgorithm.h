#ifndef CLOUDTYPEPRINTERALGORITHM_H
#define CLOUDTYPEPRINTERALGORITHM_H

#include "AbstractAlgorithm.h"
#include "GeneralConvertor.h"
#include "SingleStar.h"

class CloudTypePrinterAlgorithm
{
public:
    CloudTypePrinterAlgorithm(std::string path):path(path)
    {

    }

    void CompleteOperation(AbstractCloudType* k)
    {
        SingleStar* star1 = (SingleStar*) k;
        GeneralConvertor::FinalizeDeal(path, star1->returnNumberOfRuns(),star1->returnNumberOfBranches());
    }

protected:

private:
    std::string path;
};

#endif // CLOUDTYPEPRINTERALGORITHM_H

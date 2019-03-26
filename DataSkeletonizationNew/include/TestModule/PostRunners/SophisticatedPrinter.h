#ifndef SOPHISTICATEDPRINTER_H
#define SOPHISTICATEDPRINTER_H

#include "PostRunInterface.h"
#include "GeneralConvertor.h"
#include "SingleStar.h"

class SophisticatedPrinter : public PostRunInterface
{
public:
// Star8Collection
    SophisticatedPrinter(std::string folder): folder(folder) {}
    void run(MyGraphType & G, std::list<Point> & cloud, AbstractCloudType* gen, int RunNumber, std::string AlgorithmName) override
    {
        SingleStar* type = (SingleStar*) gen;
        int qqq = type ->returnNumberOfBranches();
        std::string subfolder = "Star" + std::to_string(qqq) + "Collection/";
        std::string finalpath = folder + subfolder + "File" + std::to_string(RunNumber) + ".xyz";
        GeneralConvertor::CloudToXYZ(cloud,finalpath,RunNumber);
       // GeneralConvertor::GraphToVtk(finalpath,G);
    }

protected:

private:
    std::string folder;
};

#endif // SOPHISTICATEDPRINTER_H

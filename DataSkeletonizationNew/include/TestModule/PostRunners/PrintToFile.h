#ifndef PRINTTOFILE_H
#define PRINTTOFILE_H

#include "PostRunInterface.h"
#include "GeneralConvertor.h"

class PrintToFile : public PostRunInterface
{
    public:
        PrintToFile(std::string folder):folder(folder) {}

        void run(MyGraphType & G, std::list<Point> & cloud, AbstractCloudType* gen, int RunNumber, std::string AlgorithmName) override
        {
        std::string finalpath = folder + gen->returnName() + AlgorithmName + std::to_string(RunNumber) + ".vtk";
        GeneralConvertor::GraphToVtk(finalpath,G);
        }

    protected:

    private:
    std::string folder;
};

#endif // PRINTTOFILE_H

#ifndef REALCLOUDCOLLECTION_H
#define REALCLOUDCOLLECTION_H

#include <boost/filesystem.hpp>
#include "AbstractCloudType.h"
#include "generatable.h"

class RealCloudCollection : public AbstractCloudType
{
public:
    RealCloudCollection(std::string nameOfInstance, std::string pathName);
    ~RealCloudCollection();
    bool IsGraphCorrect(MyGraphType & G, int iterationNumber);
    bool DoesGraphHaveCorrectForm(MyGraphType & G, int iterationNumber);
    void GenerateCloud(std::list<Point> & p, int iterationNumber);
    void SetCorrectnessOfGraphs(std::string filepath);



protected:

private:
    std::vector<boost::filesystem::directory_entry> directories;
    std::vector<generatable*> correctGraphType;
};

#endif // REALCLOUDCOLLECTION_H

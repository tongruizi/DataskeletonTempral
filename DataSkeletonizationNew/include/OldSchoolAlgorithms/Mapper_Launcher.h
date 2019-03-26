#ifndef MAPPER_LAUNCHER_H
#define MAPPER_LAUNCHER_H

#include "Graph.h"
#include "Definitions.h"
#include "Mapper.h"
#include "Mapper_Parameters.h"
#include "AbstractAlgorithm.h"

class Mapper_Launcher : public AbstractAlgorithm
{
public:
    Mapper_Launcher(Mapper_Parameters & k):mp(k),AbstractAlgorithm("Mapper")
    {}
    void Run(std::list<Point> & cloud, MyGraphType & G)
    {
        std::vector<Point> v{ std::begin(cloud), std::end(cloud) };
        Mapper(v,mp,G);
    }

protected:

private:
    Mapper_Parameters mp;

};

#endif // MAPPER_LAUNCHER_H

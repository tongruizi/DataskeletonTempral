#ifndef GRAPHGENERATION_H
#define GRAPHGENERATION_H

#include <Graph.h>
#include <Definitions.h>

class GraphGeneration
{
    public:
        GraphGeneration();
        static void RandomGraph1(int n, double minangle, double scale, MyGraphType & G);
        static void RandomGraph2(int n, int n2, double minangle, double scale,  MyGraphType & G);
        static vertex_descriptor TemplateFunction(Point extra, double xshift, double yshift,
        double zshift, MyGraphType & G, int n, double minangle, double scale);

    protected:

    private:
};

#endif // GRAPHGENERATION_H

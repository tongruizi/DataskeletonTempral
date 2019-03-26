#ifndef GRAPHSIMPLIFICATION_H
#define GRAPHSIMPLIFICATION_H

#include "Graph.h"

class GraphSimplification
{
    public:
        GraphSimplification();
        static void SimplifyGraph(MyGraphType G, MyGraphType & S);
        static bool CorrectNumberOfEnds(int k, MyGraphType & G);
        static bool RecognizeStraGraph(MyGraphType & G, int n);
        static bool RecognizeDoubleStarGraph(MyGraphType & G, int n, int m);


    protected:

    private:
};

#endif // GRAPHSIMPLIFICATION_H

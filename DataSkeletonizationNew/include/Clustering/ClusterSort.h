#ifndef CLUSTERSORT_H
#define CLUSTERSORT_H
#include "Graph.h"
#include "ClusterElement.h"

// std::pair<vertex_descriptor,double>

template <class T>


class ClusterSort
{
    public:
        ClusterSort();
        bool operator() (const ClusterElement<T>& struct1, const ClusterElement<T>& struct2);

    protected:

    private:
};

#endif // CLUSTERSORT_H

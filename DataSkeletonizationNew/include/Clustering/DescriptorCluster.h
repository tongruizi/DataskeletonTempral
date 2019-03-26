#ifndef DESCRIPTORCLUSTER_H
#define DESCRIPTORCLUSTER_H

#include "Graph.h"
#include "ClusterElement.h"

class DescriptorCluster
{
    std::pair<vertex_descriptor, double> theElement;
    public:
        DescriptorCluster(std::pair<vertex_descriptor, double> element);
        double returnValue() const;
        vertex_descriptor returnElement();
    protected:

    private:
};

#endif // DESCRIPTORCLUSTER_H

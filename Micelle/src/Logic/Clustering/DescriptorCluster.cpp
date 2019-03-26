#include "DescriptorCluster.h"
#include "Graph.h"

DescriptorCluster::DescriptorCluster(std::pair<vertex_descriptor, double> paird) :
theElement(paird)
{


}

double DescriptorCluster::returnValue() const
{
return theElement.second;
}

vertex_descriptor DescriptorCluster::returnElement()
{
return theElement.first;
}

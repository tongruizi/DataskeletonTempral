#ifndef CLUSTERINGMETHODS_H
#define CLUSTERINGMETHODS_H
#include "ClusterSort.h"

template <class T>
class ClusteringMethods
{//     int indexx = ClusteringMethods::DownwardClusterization(examElements,scale,settings);

    public:
        ClusteringMethods();
        double computeAverage(std::vector<ClusterElement<T>> & elements);
        int DownwardClusterization(std::vector<ClusterElement<T>> & elements, double coefficent, std::string preference);
        double computeDeviation(std::vector<ClusterElement<T>> & elements);
        void debugSet(std::vector<ClusterElement<T>> & elements, std::string out);

    protected:

    private:
};

#endif // CLUSTERINGMETHODS_H

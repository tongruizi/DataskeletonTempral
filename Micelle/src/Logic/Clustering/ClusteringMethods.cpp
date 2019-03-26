#include "ClusteringMethods.h"
template <class T>
ClusteringMethods<T>::ClusteringMethods()
{
    //ctor
}

template <class T>
double ClusteringMethods<T>::computeAverage(std::vector<ClusterElement<T>> & elements)
{
    double sum = 0;
    double number = 0;

    for (auto it = elements.begin(); it != elements.end(); it++)
    {
        sum = sum + (*it).returnValue();
        number = number + 1;
    }

    if (number != 0)
    {
        return sum/number;
    }
    else
    {
        return 0;
    }
}

template <class T>
double ClusteringMethods<T>::computeDeviation(std::vector<ClusterElement<T>> & elements)
{
    double number = 0;
    double sum = 0;
    double avg = ClusteringMethods<T>::computeAverage(elements);
    for (auto it = elements.begin(); it != elements.end(); it++)
    {
        double vv = (*it).returnValue();
        sum = sum + (avg-vv)*(avg-vv);
        number = number + 1;
    }
    return sqrt(sum/number);
}

template <class T>
void ClusteringMethods<T>::debugSet(std::vector<ClusterElement<T>> & elements, std::string out)
{
    std::ofstream of(out);
    of << "Average of dataset: " <<  ClusteringMethods<T>::computeAverage(elements) << std::endl;
    of << "Deviation of the dataset: " << ClusteringMethods<T>::computeDeviation(elements) << std::endl;
    for (int i = 0; i < elements.size(); i++)
    {
        of << "i: " << i << " value: " << elements[i].returnValue() << std::endl;
    }
    of.close();
}


template <class T>
int ClusteringMethods<T>::DownwardClusterization(std::vector<ClusterElement<T>> & elements, double coefficent, std::string preference)
{
//   std::cout << "Size of the container:" << elements.size() << std::endl;
    double ClusterValue = 0;
    std::sort(elements.rbegin(),elements.rend(),ClusterSort<T>());
    if (preference == "pure")
    {
        ClusterValue = coefficent;
    }
    else if (preference == "average")
    {
        double average = ClusteringMethods::computeAverage(elements);
        ClusterValue = average * coefficent;
    }
    else if (preference == "sd")
    {
        ClusterValue = coefficent * computeDeviation(elements);
    }
    int j = 0;
    for (j = elements.size()-2; j > -1; j--)
    {
        double difference = elements[j].returnValue()- elements[j+1].returnValue();
        if (difference > ClusterValue)
        {
            break;
        }
    }
    return j;
}

template class ClusteringMethods<vertex_descriptor>;
template class ClusteringMethods<edge_descriptor>;



// int j;
//    for (j = edges.size()-2; j > -2; j--)
//    {
//        if (j == -1)
//        {
//            break;
//        }
//
//        double difference = edges[j].second - edges[j+1].second;
//        if (difference > param)
//        {
//            break;
//        }
//    }
//    if (j != -1)
//    {
//        for (int k = 0; k <= j; k++)
//        {
//        boost::remove_edge(edges[k].first,mst);
//        }
//    }

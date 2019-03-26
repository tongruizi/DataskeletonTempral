#include "ClusterSort.h"

template <class T>

ClusterSort<T>::ClusterSort()
{
    //ctor
}
//bool SimplexSort::operator() (Chandler v, Chandler u)

template <class T>

bool ClusterSort<T>::operator() (const ClusterElement<T>& struct1, const ClusterElement<T>& struct2)
{
        return (struct1.returnValue() < struct2.returnValue());
}
template class ClusterSort<vertex_descriptor>;
template class ClusterSort<edge_descriptor>;

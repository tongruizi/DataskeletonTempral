#include "ClusterElement.h"
#include "Graph.h"
//
// T theElemet;
//    double value;
template <class T>


ClusterElement<T>::ClusterElement(T object, double valued):
theElement(object), value(valued)
{
    //ctor
}
template <class T>

double ClusterElement<T>::returnValue() const
{
return value;
}

template <class T>


T ClusterElement<T>::returnObject() const
{
return theElement;
}

template <class T>


bool ClusterElement<T>::operator < (const ClusterElement<T> & k) const
{
this -> returnValue() < k.returnValue();
}

template class ClusterElement<vertex_descriptor>;
template class ClusterElement<edge_descriptor>;


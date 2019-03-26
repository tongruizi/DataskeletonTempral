#ifndef CLUSTERELEMENT_H
#define CLUSTERELEMENT_H


template <class T>
class ClusterElement
{
    T theElement;
    double value;

        public:
        ClusterElement(T element, double value);
        double returnValue() const;
        T returnObject() const;
        bool operator < (const ClusterElement & k) const;
    protected:

    private:
};

#endif // CLUSTERELEMENT_H

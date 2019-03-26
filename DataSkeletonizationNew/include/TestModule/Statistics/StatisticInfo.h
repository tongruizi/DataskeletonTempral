#ifndef STATISTICINFO_H
#define STATISTICINFO_H

#include <string>
#include <limits>
#include "AbstractMeasurer.h"
//! This will be template file, so please code everything into header and completely DESTROY the implementation file

//! unit here can be either double, size_t, int or what ever. WE DONT CARE XD (but we still want stuff to work smoothly)

template<typename unit>
class StatisticInfo
{
public:
    StatisticInfo():number(0),sum(0)
    {
        maxValue=std::numeric_limits<unit>::min();
        minValue=std::numeric_limits<unit>::max();
    }
    //! wrong, we are not running anything here, this is only statistic file!
    //static void run( MyGraphType & G, std::list<Point> & points)=0;
    //! We are not selecting anything...
    //static void seleAlgorithm(std::string type);
    //! INSTEAD:
    unit returnSum()
    {
        return sum;
    }
    double returnAvg()
    {
        return (double)sum / (double) number;
    }
    unit returnMax()
    {
        return maxValue;
    }
    unit returnMin()
    {
        return minValue;
    }
    void reset()
    {
        number=0;
        sum=0;
        maxValue=std::numeric_limits<unit>::min();
        minValue=std::numeric_limits<unit>::max();
    }
    std::string returnInfo()
    {
        //! This will be implemented later
        std::string info = "";
        return info;
    }
    void add(unit t)
    {
        number = number + 1;
        sum = sum + t;
        maxValue = std::max(t, maxValue);
        minValue = std::min(t, minValue);

    }
    double returnNumber()
    {
        return number;
    }



protected:

private:

//! This is number of measurments: For now its always double, so we dont get stupid errors when dividing stuff
    double number;
//! sum is unit
    unit sum;
    unit maxValue;
    unit minValue;

};





#endif // STATISTICINFO_H

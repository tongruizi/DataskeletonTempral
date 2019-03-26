#ifndef POINTSEGMENTDISTANCE_H
#define POINTSEGMENTDISTANCE_H

#include <mlpack/core.hpp>


class PointSegmentDistance
{
public:
    PointSegmentDistance()
    {

    }
    double Calculate(arma::mat & segment, arma::vec & p)
    {
        double t = arma::dot(p-segment.col(0),segment.col(1)-segment.col(0)) / (arma::dot(segment.col(1)-segment.col(0),segment.col(1)-segment.col(0)));
        double u = std::min(std::max(t,(double) 0), (double) 1);
        return arma::norm(segment.col(0)+u*(segment.col(1)-segment.col(0))-p);
        //  return 0;

    }
protected:

private:
};

#endif // POINTSEGMENTDISTANCE_H

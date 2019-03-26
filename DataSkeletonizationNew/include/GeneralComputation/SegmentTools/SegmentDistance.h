#ifndef SEGMENTDISTANCE_H
#define SEGMENTDISTANCE_H

#include "mlpack/core/metrics/lmetric.hpp"
#include "Definitions.h"

class SegmentDistance
{
public:
    SegmentDistance() {}
    template<typename VecTypeA, typename VecTypeB>
    static typename VecTypeA::elem_type Evaluate(const VecTypeA& a,
            const VecTypeB& b)
    {
      //  std::cout << "S:" << Point(a(0),a(1),a(2)) << " - " << Point(a(3),a(4),a(5)) <<
      //  "ctS: " << Point(b(0),b(1),b(2)) << " - " << Point(b(3),b(4),b(5)) << std::endl;
        double d = sqrt(CGAL::squared_distance(Segment(Point(a(0),a(1),a(2)),Point(a(3),a(4),a(5))),Segment(Point(b(0),b(1),b(2)),Point(b(3),b(4),b(5)))));
     //   std::cout << "Resulting into: " << d << std::endl;
        return d;

    }


protected:

private:
};

#endif // SEGMENTDISTANCE_H

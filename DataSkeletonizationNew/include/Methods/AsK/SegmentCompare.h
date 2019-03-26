#ifndef SEGMENTCOMPARE_H
#define SEGMENTCOMPARE_H

#include "Definitions.h"

struct SegmentCompare
{
    bool operator()(const Segment& a, const Segment& b) const
    {
        if (a.source() == b.source())
        {
            return a.target() < b.target();
        }
        return a.source() < b.source();
    }
};
#endif // SEGMENTCOMPARE_H

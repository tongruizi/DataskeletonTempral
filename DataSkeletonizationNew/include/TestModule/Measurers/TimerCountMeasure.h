#ifndef TIMERCOUNTMEASURE_H
#define TIMERCOUNTMEASURE_H

#include "Definitions.h"
#include "ManualMeasure.h"
#include "AbstractCloudType.h"

class TimerCountMeasure
{
public:
    TimerCountMeasure(std::vector<std::pair<std::string,std::string>> & timers):
        timers(timers),sums(timers.size())
    {
        std::fill(sums.begin(), sum.end(), 0);

    }

    std::string ReturnTimers(int q, long numberOfRuns)
    {
        std::string s = "";

        for(int i = 0; i < timers.size(); i++)
        {
            long midr = mlpack::Timer::Get(timers[i].second + q).count();
            long avg = midr / numberOfRuns;
            if (i < timers.size() - 1)
            {
                s = s + timers[i].first + avg + " + ";
            }
            else
            {
                s = s + timers[i].first + avg;
            }
        }

        return s;
    }

    TimerCountMeasure* Clone() override
    {
    return new TimeCountMeasure(*this);
    }

protected:

private:
    std::vector<std::pair<std::string,std::string>> timers;
};

#endif // TIMERCOUNTMEASURE_H

#ifndef DOUBLESTAR_H
#define DOUBLESTAR_H

#include "CloudGenerator.h"
#include "GraphGeneration.h"
#include "generatable.h"
#include <set>

class DoubleStar:public generatable
{
public:
    //! Implement this, as in single star. Note that we have an extra parameter...
    DoubleStar(double angle, int number_of_branches,int number_of_branches2, int number_of_cloudpoints, double epsilon, double scale,int number_of_runs,
    std::string name = "DoubleStar"):
        generatable(number_of_cloudpoints,epsilon,number_of_runs,name),angle(angle), number_of_branches(number_of_branches),
         number_of_branches2(number_of_branches2),scale(scale)
    {

    }

    void GenerateGraph(MyGraphType & G)
    {
        GraphGeneration::RandomGraph2(number_of_branches, number_of_branches2, angle, scale, G);
    }

    bool IsGraphCorrect(MyGraphType & G, int iterationnumber)
    {
        std::set<int> setInt;
        int n2 = number_of_branches + 1;
        int m2 = number_of_branches2 + 1;
        setInt.insert(n2);
        setInt.insert(m2);
        int onev;
        auto vpair = boost::vertices(G);
        for (auto it = vpair.first; it != vpair.second; it++)
        {
            int tmpD = boost::degree(*it, G);
            if (tmpD > 2)
            {
                auto bt = setInt.find(tmpD);
                if (bt == setInt.end())
                {
                    std::cout << tmpD << std::endl;
                    return false;
                }
                else
                {
                    setInt.erase(tmpD);
                }
            }
            else if (tmpD == 1)
            {
                onev++;
            }
        }
        return (setInt.size() == 0) && (onev == number_of_branches + number_of_branches2);
    }

    bool DoesGraphHaveCorrectForm(MyGraphType & G, int iterationnumber)
    {
        return this->CorrectNumberOfEndPoints(G,number_of_branches+number_of_branches2);
    }

private:
//! Variables here:
    double angle;
    int number_of_branches;
    int number_of_branches2;
    double scale;

};

#endif // DOUBLESTAR_H

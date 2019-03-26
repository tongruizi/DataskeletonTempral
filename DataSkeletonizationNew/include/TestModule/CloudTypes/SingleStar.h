#ifndef SINGLESTAR_H
#define SINGLESTAR_H

#include "CloudGenerator.h"
#include "GraphGeneration.h"
#include "generatable.h"

class SingleStar:public generatable
{
public:
    SingleStar(double angle, int number_of_branches, int number_of_cloudpoints, double epsilon, double scale,int number_of_runs,
               std::string name = "SingleStar"):
        generatable(number_of_cloudpoints,epsilon,number_of_runs,name),angle(angle), number_of_branches(number_of_branches),scale(scale)
    {

    }
    void GenerateGraph(MyGraphType & G)
    {
        GraphGeneration::RandomGraph1(number_of_branches,angle,scale,G);
    }

    bool IsGraphCorrect(MyGraphType & G, int iterationnumber)
    {
        int n = number_of_branches;
        auto vpair = boost::vertices(G);
        int degree1 = 0;
        int starn = 0;
        for (auto it = vpair.first; it != vpair.second; it++)
        {
            int tmpD = boost::degree(*it, G);
            if (tmpD == 1)
            {
                degree1++;
            }
            else if (tmpD == n)
            {
                starn++;
            }
            else if ((tmpD != 2) && (tmpD != n) && (tmpD != 1))
            {
                // std::cout << "We exited the method because we encountered degree: " << tmpD << std::endl;
                return false;
            }
        }
        //   std::cout << "Count of the star vertices " << starn << " | " << "count of degree 1 vertices: " << degree1 << std::endl;
        if (n == 2)
        {
            return (degree1 == number_of_branches);

        }
        else
        {
            return (starn == 1) && (degree1 == number_of_branches);
        }
    }

    bool DoesGraphHaveCorrectForm(MyGraphType & G,int iterationnumber)
    {
        bool test = this->CorrectNumberOfEndPoints(G,number_of_branches);
        if (test == false)
        {
            std::cout << "INCORRECT NUMBER OF ENDPOINTS!!!!!!!" << std::endl;
        }
        return test;
    }
    int returnNumberOfBranches()
    {
        return number_of_branches;
    }



private:
    double angle;
    int number_of_branches;
    double scale;
};

#endif // SINGLESTAR_H

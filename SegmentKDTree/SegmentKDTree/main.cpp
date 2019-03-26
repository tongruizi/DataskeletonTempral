#include <iostream>
#include <mlpack/core.hpp>
#include "PointSegmentDistance.h"
#include <stdio.h>
#include <stdlib.h>
#include <mlpack/core/tree/binary_space_tree/typedef.hpp>


using namespace std;

void TestOne()
{
    arma::mat k(3,1);
    std::cout << k << std::endl;

}

void TestTwo()
{

    PointSegmentDistance dm;
    arma::vec p = {-3,1,0};
    arma::mat k(3,2);
    arma::vec v = {-1,1,0};
    arma::vec u = {1,1,0};
    k.col(0) = v;
    k.col(1) = u;
//k.replace(k.col(0), {-1,0,0});
//k.replace(k.col(1), {1,0,0});
    double result = dm.Calculate(k,p);
//std:cout << k << std::endl;
    std::cout << "Result: " << result << std::endl;
}

void TestDebugger()
{
    arma::mat k(2,25);
    int ct = 0;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            arma::vec v = {i,j};
            k.col(ct) = v;
            ct++;
        }
    }
    std::vector<size_t> oldFromNew;
    mlpack::tree::KDTree<mlpack::metric::EuclideanDistance,mlpack::tree::EmptyStatistic,arma::mat> w(k,oldFromNew);

    std::cout << k << std::endl;
    std::cout << "Succeful" << std::endl;

}

int main()
{
TestDebugger();
}

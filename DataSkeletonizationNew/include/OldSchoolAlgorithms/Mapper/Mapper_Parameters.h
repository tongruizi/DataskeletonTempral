#pragma once

#include <string>


class Mapper_Parameters
{
    public:

    int num_intervals;
    double mcsf, overlap_ratio, sigma;
    std::string filter_function;

    Mapper_Parameters ( int n_i, double o_r, std::string const& f_f, double s, double m );

    Mapper_Parameters();
    ~Mapper_Parameters();
};

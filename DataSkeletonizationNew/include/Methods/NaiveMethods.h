#ifndef NAIVEMETHODS_H
#define NAIVEMETHODS_H

#include <mlpack/core.hpp>
#include <mlpack/prereqs.hpp>


class NaiveMethods
{
public:
    NaiveMethods();
    template<typename f>
    static void ComputeDensityNaiveRow(std::vector<double> & results, arma::mat & data)
    {
        results.resize(data.n_rows);
        for (int i = 0; i < data.n_rows; i++)
        {
            for (int j = 0; j < data.n_rows; j++)
            {
                results[i] = results[i] + f::Evaluate(mlpack::metric::EuclideanDistance::Evaluate(data.row(i), data.row(j)));
            }
        }
    }
     template<typename f>
    static void ComputeDensityNaiveCol(std::vector<double> & results, arma::mat & data)
    {
        results.resize(data.n_cols);
        for (int i = 0; i < data.n_cols; i++)
        {
            for (int j = 0; j < data.n_cols; j++)
            {
                results[i] = results[i] + f::Evaluate(mlpack::metric::EuclideanDistance::Evaluate(data.col(i), data.col(j)));
            }
        }
    }
protected:

private:
};

#endif // NAIVEMETHODS_H

#ifndef POINTEXTRACTOR_H
#define POINTEXTRACTOR_H

#include <mlpack/core.hpp>

class PointExtractor
{
public:
    PointExtractor() {}
    //! Please use this only if Segments is non empty
    static void ExtractPoints(arma::mat & dataset, std::vector<arma::mat> & Segments)
    {
        if (Segments.size() == 0)
        {
            return;
        }
        int sgSize = Segments.size();
        int dim = Segments[0].n_rows;
        dataset.set_size(dim,2*sgSize);
        for (int i = 0; i < Segments.size(); i++)
        {
            dataset.col(i) = Segments[i].col(0);
            dataset.col(i + sgSize) = Segments[i].col(1);
        }
    }
    static void ExtractPoints(arma::mat & dataset, std::vector<arma::mat> & Segments, std::vector<size_t> & DataToSegment)
    {
        if (Segments.size() == 0)
        {
            return;
        }
        int sgSize = Segments.size();
        int dim = Segments[0].n_rows;
        dataset.set_size(dim,2*sgSize);
        DataToSegment.set_size(2*sgSize);
        for (int i = 0; i < Segments.size(); i++)
        {
            dataset.col(i) = Segments[i].col(0);
            DataToSegment[i] = i;
            dataset.col(i + sgSize) = Segments[i].col(1);
            DataToSegment[i + sgSize] = i;
        }
    }
    //  virtual ~PointExtractor() {}

protected:

private:
};

#endif // POINTEXTRACTOR_H

#ifndef GENERALCONVERTOR_H
#define GENERALCONVERTOR_H

#include <string>
#include <mlpack/core.hpp>
#include <vector>
#include "Graph.h"
#include "Definitions.h"

//! Matinfo to file:


class GeneralConvertor
{
public:
    GeneralConvertor();
    static void XYZtoMAT(std::string fl, arma::mat & data);
    static void XYZtoPoint(std::string fl, std::list<Point> & points);
    static void MatInfoToFile(std::string out, arma::mat & data, std::vector<double> & scalar);
    static void VectorToFile(std::string out, std::vector<int> & s);
    static void MSTToVTK(arma::mat & originalData, arma::mat & edges, std::string out);
    static void GraphToVtk(std::string path, MyGraphType &G);
    static void ListToMat(std::list<Point> & lp, arma::mat & mtx);
    static void ListToMatTransposed(std::list<Point> & points, arma::mat & data);
    static double ArmaMatToGraph(MyGraphType & G, arma::mat & edges, arma::mat & originalData);
    static void MatToMyGraphType(arma::mat & originalData, arma::mat & edges, MyGraphType & G);
    static void pathPrintToVtkPointlist(std::list<std::list<Point>> & paths, std::string directory);
    static void ClusteringInfoToFile(std::list<Point> & cloud, arma::Mat<size_t> & theNeighbors, std::string k,int graphsize);
    static void GraphToPaths(MyGraphType & G, std::vector<Segment> & segments);
    static void SegmentsToMat(std::vector<Segment> & segments, arma::mat & segmentmat);




    //! Method to print data from measurers directly into latex table.
    static void DataToLatex(std::vector<std::vector<std::vector<std::string>>> & measurers,
                            std::vector<std::vector<std::string>> & timeMeasures,
                            std::vector<std::string> & GraphNames, std::vector<std::string> & AlgorithmNames,
                            std::vector<std::string> & MeasurerNames, std::string filename);
    static void DataToCSV(std::vector<std::vector<std::vector<std::string>>> & measurers,
                                     std::vector<std::vector<std::string>> & timeMeasures,
                                     std::vector<std::string> & GraphNames, std::vector<std::string> & AlgorithmNames,
                                     std::vector<std::string> & MeasureNames, std::string filename);
    static void RetriveGraphInformation(std::string filename, std::vector<int> & iterationanumber,
   std::vector<std::string> & graphType, std::vector<std::vector<int>> & parameters);

    static void CloudToXYZ(std::list<Point> & p, std::string fp, int i);
    static void FinalizeDeal(std::string fp, int sz, int k);
    static void StraighteningDebugPrint(std::string filename,  std::vector<std::vector<Point>> & debuglist);
    // :ArmaMatToGraph
    static void VectorToMatTransposed(std::vector<Point> & points, arma::mat & data);





protected:

private:
};

#endif // GENERALCONVERTOR_H

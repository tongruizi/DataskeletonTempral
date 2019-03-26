#include "GeneralConvertor.h"
#include <fstream>
#include "Graph.h"

GeneralConvertor::GeneralConvertor()
{
    //ctor
}

void GeneralConvertor::XYZtoMAT(std::string fl, arma::mat & data)
{
    // We make sure that the data-matrix is empy
    data.reset();
    std::ifstream infile(fl.c_str());
    double x;
    double y;
    double z;
    std::string w;
    int k;
    infile >> k;
    data.set_size(k,3);
    getline(infile,w);
    getline(infile,w);
    //  std::vector<std::vector<double>> toArma;
    int i = 0;
    while (infile >> w >> x >> y >> z)
    {
        data(i,0) = x;
        data(i,1) = y;
        data(i,2) = z;
        i++;
    }
}

void GeneralConvertor::XYZtoPoint(std::string fl, std::list<Point> & points)
{
    // We make sure that the data-matrix is empy
    std::ifstream infile(fl.c_str());
    double x;
    double y;
    double z;
    std::string w;
    int k;
    infile >> k;
    getline(infile,w);
    getline(infile,w);
    //  std::vector<std::vector<double>> toArma;
    int i = 0;
    while (infile >> w >> x >> y >> z)
    {
        points.push_back(Point(x,y,z));
        i++;
    }

}



void GeneralConvertor::MatInfoToFile(std::string out, arma::mat & data, std::vector<double> & scalar)
{
    std::ofstream mystream;
    mystream.open(out);
    mystream << "x coord, y coord, z coord, scalar" << std::endl;
    for (int i = 0; i < data.n_rows; i++)
    {
        double v1 = data(i,0);
        double v2 = data(i,1);
        double v3 = data(i,2);
        mystream << v1;
        mystream << ", ";
        mystream << v2;
        mystream << ", ";
        mystream << v3;
        mystream << ", ";
        mystream << scalar[i];
        mystream << std::endl;
        // ", " << v2 ", " << v3 << ", " << scalar[i] << std::endl;
    }
    mystream.close();
}

void GeneralConvertor::ClusteringInfoToFile(std::list<Point> & cloud, arma::Mat<size_t> & theNeighbors, std::string out, int graphsize)
{
    std::vector<int> shuffleFunction;
    for (int i = 0; i < graphsize; i++)
    {
        shuffleFunction.push_back(i);
    }
    std::random_shuffle(shuffleFunction.begin(),shuffleFunction.end());
    std::ofstream mystream;
    mystream.open(out);
    mystream << "x coord, y coord, z coord, scalar" << std::endl;
    size_t number = 0;
    for (auto it = cloud.begin(); it != cloud.end(); it++)
    {
        double cl = (double) shuffleFunction[theNeighbors(0,number)]/(double) graphsize;
        mystream << (*it).x() << ", " << (*it).y() << ", " << (*it).z() << ", " << cl << std::endl;
        number++;
    }

    mystream.close();

}

void GeneralConvertor::VectorToFile(std::string out, std::vector<int> & s)
{

    std::ofstream mystream;
    mystream.open(out);
    for (int i = 0; i < s.size(); i++)
    {
        mystream << i << ", " << s[i] << std::endl;
    }
    mystream.close();

}
void GeneralConvertor::pathPrintToVtkPointlist(std::list<std::list<Point>> & paths, std::string directory)
{

    int size = 0;
    for (auto iter = paths.begin(); iter != paths.end(); iter++)
    {
        size = size + (*iter).size();
    }
    std::ofstream mystream;
    mystream.open(directory);
    mystream << "# vtk DataFile Version 1.0\n";
    mystream << "3D triangulation data\n";
    mystream << "ASCII\n";

    mystream << std::endl;
    mystream << "DATASET POLYDATA\n";
    mystream << "POINTS " << size << " float\n";

    for(auto globalit = paths.begin(); globalit != paths.end(); globalit++)
    {
        std::list<Point> path = *globalit;
        for (auto pathit = path.begin(); pathit != path.end(); pathit++)
        {
            mystream << (*pathit) << std::endl;
        }
    }

    mystream << "LINES " << (size-paths.size()) << " " << (size-paths.size())*3 << std::endl;

    int location = 0;
    for(auto globalit = paths.begin(); globalit != paths.end(); globalit++)
    {
        int cursize = (*globalit).size();
        for (int i = location; i < location + cursize-1; i++)
        {
            mystream << "2 " << i << " " << i+1 << std::endl;
        }
        location = location + cursize;

    }
    mystream.close();
}


void GeneralConvertor::MatToMyGraphType(arma::mat & originalData, arma::mat & edges, MyGraphType & G)
{
    int vSize = originalData.n_cols;
    for (int i = 0; i < vSize; i++)
    {
        Graph::add_vertex(G,Point(originalData(0,i),originalData(1,i),originalData(2,i)));
    }
    int edgeSize = edges.n_cols;
    for (int i = 0; i < edgeSize; i++)
    {
        Graph::add_edge(G,edges(0,i),edges(1,i));
    }

}

void GeneralConvertor::MSTToVTK(arma::mat & originalData, arma::mat & edges, std::string out)
{
    std::ofstream mystream;
    mystream.open(out);
    mystream << "# vtk DataFile Version 1.0" << std::endl;
    mystream << "3D triangulation data" << std::endl;
    mystream << "ASCII" << std::endl;
    mystream << std::endl;
    mystream << "DATASET POLYDATA" << std::endl;
    int vSize = originalData.n_cols;
    mystream << "POINTS " << vSize << " float" << std::endl;
    for (int i = 0; i < vSize; i++)
    {
        mystream << originalData(0,i) << " " << originalData(1,i) << " " << originalData(2,i) << std::endl;
    }
    int edgeSize = edges.n_cols;
    int tripledgeSize = 3 * edgeSize;
    mystream << "LINES " << edgeSize << " " << tripledgeSize << std::endl;
    for (int i = 0; i < edgeSize; i++)
    {
        mystream << "2 " << edges(0,i) << " " << edges(1,i) << std::endl;
    }
    mystream.close();
}

void GeneralConvertor::GraphToVtk(std::string path, MyGraphType &G)
{
    std::ofstream mystream;
    mystream.open(path);
    mystream << "# vtk DataFile Version 1.0\n";
    mystream << "3D triangulation data\n";
    mystream << "ASCII\n";
    mystream << std::endl;
    mystream << "DATASET POLYDATA\n";

    mystream << "POINTS " << num_vertices(G) << " float\n";
    for(int i=0; i<num_vertices(G); i++)
    {
        mystream << G[i].p << std::endl;
    }

    mystream << "LINES " << (num_edges(G)) << " " << (num_edges(G))*3 << std::endl;

    auto epair = edges(G);
    for(auto iter=epair.first; iter!=epair.second; iter++)
    {
        mystream  << "2 " << source(*iter, G) << " " << target(*iter, G) <<std::endl;
    }
    mystream.close();
}

void GeneralConvertor::GraphToPaths(MyGraphType & G, std::vector<Segment> & segments)
{
    auto epair = boost::edges(G);
    for (auto it = epair.first ; it != epair.second ; it++)
    {
        segments.push_back(Segment(G[source(*it,G)].p, G[target(*it,G)].p));
    }


}



void GeneralConvertor::ListToMat(std::list<Point> & points, arma::mat & data)
{
    int k = points.size();
    data.set_size(k,3);
    int j = 0;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        data(j,0) = it->x();
        data(j,1) = it->y();
        data(j,2) = it->z();
        j++;
    }
}
void GeneralConvertor::ListToMatTransposed(std::list<Point> & points, arma::mat & data)
{
    int k = points.size();
    data.set_size(3,k);
    int j = 0;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        data(0,j) = it->x();
        data(1,j) = it->y();
        data(2,j) = it->z();
        j++;
    }

}
double GeneralConvertor::ArmaMatToGraph(MyGraphType & G, arma::mat & edges, arma::mat & originalData)
{
    int vSize = originalData.n_cols;
    for (int i = 0; i < vSize; i++)
    {
        Point p(originalData(0,i),originalData(1,i),originalData(2,i));
        Graph::add_vertex(G,p);
    }

    int edgeSize = edges.n_cols;
    double sum = 0;

    for (int i = 0; i < edgeSize; i++)
    {
        int one = (int) edges(0,i);
        int two = (int) edges(1,i);
        Graph::add_edge(G,one,two);
        sum = sum + edges(2,i);
    }
    return sum / (double) edgeSize;

}

void GeneralConvertor::RetriveGraphInformation(std::string filename, std::vector<int> & iterationanumber,
        std::vector<std::string> & graphType, std::vector<std::vector<int>> & parameters)
{
    std::cout << "Opened file: " << filename << std::endl;
    int NumberOfLines;
    std::ifstream infile(filename);
    std::string line;
    std::string graphnamme;
    int number;
//   infile >> NumberOfLines;
//   ititerationanumber[j];erationanumber.resize(NumberOfLines);
//   graphType.resize(NumberOfLines);
    //  parameters.resize(NumberOfLines);
    int j = 0;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> number;
        iss >> graphnamme;
        graphType[j] = graphnamme;
        iterationanumber[j] = number;
        int i;
        while(iss >> i)
        {
            parameters[j].push_back(i);
        }
        j++;
    }
}

void GeneralConvertor::DataToLatex(std::vector<std::vector<std::vector<std::string>>> & measurers,
                                   std::vector<std::vector<std::string>> & timeMeasures,
                                   std::vector<std::string> & GraphNames, std::vector<std::string> & AlgorithmNames,
                                   std::vector<std::string> & MeasureNames, std::string filename)
{
    std::ofstream mystream;
    mystream.open(filename);
//    mystream << '\';
    mystream << "begin{tabular}{";
    int sum = 3 + MeasureNames.size();
    for (int i = 0; i < sum ; i++)
    {
        if (i < sum - 1)
        {
            mystream << "c|";
        }
        else
        {
            mystream << "c";
        }
    }
    mystream << "}" << std::endl;
    //  mystream << '\';
    mystream << "toprule" << std::endl;

    // graph \, & \, algorithm \, & \, endpoints \, & \, homeo type \, & \, time, ms \, & \, error \,
    mystream << "graph \\, & \\, algorithm \\, & \\, time \\, ";
    for (int i = 0; i < MeasureNames.size(); i++)
    {
        mystream << "& \\, " << MeasureNames[i] << " \\,";
    }
    mystream << std::endl;
    mystream << "\\\\" << std::endl;
    for (int i = 0; i < GraphNames.size(); i++)
    {
        mystream << "\\midrule" << std::endl;
        for (int j = 0; j < AlgorithmNames.size(); j++)
        {
            // 3-star & Mapper & 75\% & 75\% & 611 & 18.9\% \\ %& 18.2 \\

            mystream << GraphNames[i] ;
            mystream << " & " ;
            mystream << AlgorithmNames[j] ;
            mystream << " & ";
            mystream << timeMeasures[i][j];
            mystream << "%";
            for (int k = 0; k < MeasureNames.size(); k++)
            {
                mystream << " & ";
                mystream << measurers[i][j][k];
                mystream << "%";

            }
            mystream << "\\\\" << std::endl;

        }

    }
    mystream << "\\bottomrule" << std::endl;
    mystream << "\\end{tabular}";
    mystream.close();

}

void GeneralConvertor::DataToCSV(std::vector<std::vector<std::vector<std::string>>> & measurers,
                                 std::vector<std::vector<std::string>> & timeMeasures,
                                 std::vector<std::string> & GraphNames, std::vector<std::string> & AlgorithmNames,
                                 std::vector<std::string> & MeasureNames, std::string filename)
{
    std::ofstream mystream;
    mystream.open(filename);
    mystream << "Graph, Algorithm, Time (s), ";
    for (int i = 0; i < MeasureNames.size(); i++)
    {
        mystream << MeasureNames[i];

        if (i < MeasureNames.size() - 1)
        {
            mystream << ", ";
        }
    }
    mystream << std::endl;
    for (int i = 0; i < GraphNames.size(); i++)
    {
        for (int j = 0; j < AlgorithmNames.size(); j++)
        {
            // 3-star & Mapper & 75\% & 75\% & 611 & 18.9\% \\ %& 18.2 \\

            mystream << GraphNames[i] ;
            mystream << ", " ;
            mystream << AlgorithmNames[j] ;
            mystream << ", ";
            mystream << timeMeasures[i][j];
            for (int k = 0; k < MeasureNames.size(); k++)
            {
                mystream << ", ";
                mystream << measurers[i][j][k];
            }
            mystream << std::endl;
        }
    }
    mystream.close();
}

void GeneralConvertor::CloudToXYZ(std::list<Point> & p, std::string fp, int i)
{
    std::ofstream mystream;
    mystream.open(fp);
    std::string tmp1 = std::to_string(p.size());
    mystream << tmp1 << std::endl;
    mystream << std::to_string(i) << std::endl;
    int j = 0;
    for (auto it = p.begin(); it != p.end(); it++)
    {
        mystream << std::to_string(j) << " " << (*it).x() << " " << (*it).y() << " " << (*it).z() << std::endl;
        j++;
    }
    mystream.close();
}

void GeneralConvertor::FinalizeDeal(std::string fp, int sz, int k)
{
    std::ofstream mystream;
    mystream.open(fp);
    for (int i = 0; i < sz-1; i++)
    {
        mystream << i << " SingleStar " << k << std::endl;
    }
    int kam = sz-1;
    mystream << kam << " SingleStar " << k;
    mystream.close();

}


void GeneralConvertor::StraighteningDebugPrint(std::string filename,  std::vector<std::vector<Point>> & debuglist)
{

    std::ofstream mystream;
    mystream.open(filename);
    mystream << "x coord, y coord, z coord, scalar" << std::endl;

    for (int i = 0; i < debuglist.size(); i++)
    {
        int sizem = debuglist[i].size();
        for (int j = 0; j < sizem; j++)
        {
            double valuev = (double) j / (double) sizem;
            mystream << debuglist[i][j].x() << ", " <<  debuglist[i][j].y() << ", " <<  debuglist[i][j].z() << ", " << valuev << std::endl;
        }
    }
    mystream.close();
}

void GeneralConvertor::SegmentsToMat(std::vector<Segment> & segments, arma::mat & segmentmat)
{

    segmentmat.set_size(6,segments.size());
    int indx = 0;
    for (int i = 0; i < segments.size(); i++)
    {
        Point p = (segments[i]).source();
        Point q = (segments[i]).target();
        segmentmat(0,i) = p.x();
        segmentmat(1,i) = p.y();
        segmentmat(2,i) = p.z();
        segmentmat(3,i) = q.x();
        segmentmat(4,i) = q.y();
        segmentmat(5,i) = q.z();
    }
}

void GeneralConvertor::VectorToMatTransposed(std::vector<Point> & points, arma::mat & data)
{
    int k = points.size();
    data.set_size(3,k);
    int j = 0;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        data(0,j) = it->x();
        data(1,j) = it->y();
        data(2,j) = it->z();
        j++;
    }

}









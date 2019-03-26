#include "UIHandler.h"
#include "CloudGenerator.h"
#include "GraphGeneration.h"
#include "Definitions.h"
#include "Computation.h"
#include "Write.h"
#include "FlexibleComplex.h"
#include "BranchDetection.h"
#include "HashComplex.h"
#include "ComplexAlgorithms.h"
#include "AlphaReebComputation.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include "Mapper.h"
#include "Filereader.h"
#include "ClassicUIHandle.h"
#include "GraphSimplification.h"

UIHandler::UIHandler()
{
    //ctor
}

void UIHandler::setFolder(std::string newFolder)
{
    folder = newFolder;
}


void AlphaBoundaryFaceInfo(Alpha_shape_3 & as)
{
    int interior = 0;
    int exterior = 0;
    int regular = 0;
    int singular = 0;
    for (auto edge = as.finite_edges_begin(); edge != as.finite_edges_end(); edge++)
    {
        if (as.classify(*edge) == Alpha_shape_3::REGULAR)
        {
            regular++;
        }
        if (as.classify(*edge) == Alpha_shape_3::INTERIOR)
        {
            interior++;
        }
        if (as.classify(*edge) == Alpha_shape_3::EXTERIOR)
        {
            exterior++;
        }
        if (as.classify(*edge) == Alpha_shape_3::SINGULAR)
        {
            singular++;
        }
    }
    std::cout << "Interior: " << interior << std::endl;
    std::cout << "Exterior: " << exterior << std::endl;
    std::cout << "Regular: " << regular << std::endl;
    std::cout << "Singular: " << singular << std::endl;



}


void AlphaBoundaryInformation2(Alpha_shape_3 & as )
{
    std::vector<Alpha_shape_3::Cell_handle> cells;
    std::vector<Alpha_shape_3::Cell_handle> cells1;
    std::vector<Alpha_shape_3::Cell_handle> cells2;

    std::vector<Alpha_shape_3::Cell_handle> cells3;

    as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_cells(std::back_inserter(cells1), Alpha_shape_3::INTERIOR);
    as.get_alpha_shape_cells(std::back_inserter(cells2), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_cells(std::back_inserter(cells3), Alpha_shape_3::EXTERIOR);


    std::cout << "Regular cells: "<< cells.size() << std::endl;
    std::cout << "Interior cells: " << cells1.size() << std::endl;
    std::cout << "Singular cells: " << cells2.size() << std::endl;
    std::cout << "Exterior cells: " << cells3.size() << std::endl;

}

void AlphaBoundaryInformation(Alpha_shape_3 & as )
{
    std::cout << "Alpha boundary operation: " << std::endl;
    std::vector<Alpha_shape_3::Facet> facets;
    std::vector<Alpha_shape_3::Edge> edges;
    std::vector<Alpha_shape_3::Vertex_handle> vertexes;

    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertexes),Alpha_shape_3::REGULAR);

    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertexes), Alpha_shape_3::SINGULAR);

    int Euler = vertexes.size() - edges.size() + facets.size();
    std::cout << "Real Euler: " << Euler << std::endl;


}

void AlphaInformation(Alpha_shape_3 & as, Alpha_iterator opt )
{

    std::vector<Alpha_shape_3::Cell_handle> cells;
    std::vector<Alpha_shape_3::Facet> facets;
    std::vector<Alpha_shape_3::Edge> edges;
    std::vector<Alpha_shape_3::Vertex_handle> vertexes;

    as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertexes),Alpha_shape_3::REGULAR);

    as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::SINGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertexes), Alpha_shape_3::SINGULAR);

    as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::INTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::INTERIOR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::INTERIOR);
    as.get_alpha_shape_vertices(std::back_inserter(vertexes), Alpha_shape_3::INTERIOR);

    std::cout << "Information of the Alpha_Complex:\n";
    std::cout << "The alpha-complex has : " << std::endl;
    std::cout << cells.size() << " cells as tetrahedrons" << std::endl;
    std::cout << facets.size() << " triangles" << std::endl;
    std::cout << edges.size() << " edges" << std::endl;
    std::cout << vertexes.size() << " vertices" << std::endl;
    int Euler = vertexes.size() - edges.size() + facets.size() - cells.size();
    std::cout << "Euler: " << Euler << std::endl;


}

void GenerateStarCloud(std::list<Point> & out, int n, int k,double minangle, double scale, double epsilon)
{
    MyGraphType G;
    GraphGeneration::RandomGraph1(k,minangle,scale,G);
    CloudGenerator::generatePoints(n, G, epsilon, out);
}

void UIHandler::AlphaMST(std::list<Point> cloud, std::string f1, std::string f2, std::string f3, std::string f4)
{
    Alpha_shape_3 as(cloud.begin(),cloud.end());
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
    //  AlphaInformation(as,opt);
    // AlphaBoundaryInformation2(as);
    AlphaBoundaryFaceInfo(as);
    MyGraphType tree = Computation::computeMST(cloud);
    Write::GraphToVtk(f1,tree);
    Write::AlphaVTK(f2, as);
    MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as,f3);
    Write::GraphToVtk(f3,ABC);
    MyGraphType optiout;
    BranchDetection::SimplifyIt(ABC,optiout,200.0,"","pure");
    Write::GraphToVtk(f4,optiout);
    //Write::BoundaryVerticesVTK(f2, as);
}

void UIHandler::runTest(int n, int k, double minangle, double scale, std::string filename, double epsilon, std::string f2, std::string f3, std::string f4)
{
    std::list<Point> cloud;
    GenerateStarCloud(cloud, n, k, minangle, scale,epsilon);
    UIHandler::AlphaMST(cloud, folder + filename, folder + f2, folder + f3, folder + f4);
}

void UIHandler::runTestPremium(int n, int k, double minangle, double scale, std::string filename,
                               double epsilon, std::string f2, std::string fi)
{
    std::list<Point> cloud;
    GenerateStarCloud(cloud, n, k, minangle, scale,epsilon);
    MyGraphType tree = Computation::computeMST(cloud);
    Computation::FiltrationOfAlphaShapes(cloud, folder + f2);
    Write::GraphToVtk(folder + fi,tree);
}

void UIHandler::multipleTestsPremium(int n, int k, double minangle, double scale, int testN, double epsilon)
{
    for (int i = 0; i < testN; i++)
    {
        std::cout << "Round: " << i << " on the way..." << std::endl;
        //  std::string filename = "mst" + std::to_string(i) + ".vtk";
        std::string f1 = "alpha" + std::to_string(i) + "/";
        std::string fi = f1 + "mst.vtk";
        std::string filename = "";
        UIHandler::runTestPremium(n,k,minangle, scale, filename,epsilon,f1,fi);
    }

}

void UIHandler::multipleTests(int n, int k, double minangle, double scale, int testN, double epsilon)
{
    for (int i = 0; i < testN; i++)
    {
        std::cout << "Round: " << i << " on the way..." << std::endl;
        std::string filename = "mst" + std::to_string(i) + ".vtk";
        std::string f2 = "alpha" + std::to_string(i) + ".vtk";
        std::string f3 = "optimst" + std::to_string(i) + ".vtk";
        std::string f4 = "optipath" + std::to_string(i) + ".vtk";
        runTest(n,k,minangle, scale, filename,epsilon,f2,f3,f4);
    }

}

void UIHandler::AlphaSimplification(std::list<Point> cloud, std::string f1, std::string f2)
{
    Alpha_shape_3 as(cloud.begin(),cloud.end());
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
    AlphaBoundaryFaceInfo(as);
    Write::AlphaVTK(f1, as);
    AbstractComplex theComplex;
    theComplex.Convert(as);
    theComplex.SimplifyComplex();
    std::vector<Point>* points = theComplex.returnVertexIndexation();
    std::vector<Chandler> goodSimplices;
    theComplex.CollectIrreducableCells(goodSimplices);

    ST* st = theComplex.returnST();

    Write::AlphaVTKSpecial(f2,goodSimplices,*points,*st);
    //Write::BoundaryVerticesVTK(f2, as);
}

void UIHandler::AlphaDistanceComputer(std::list<Point> cloud, std::string f1, std::string f2, std::string f3, std::string f4)
{
    Alpha_shape_3 as(cloud.begin(),cloud.end());
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
    MyGraphType tree = Computation::computeMST(cloud);
    Write::GraphToVtk(f1,tree);
    Write::AlphaVTK(f2, as);
    MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as, f3);
    Write::GraphToVtk(f3,ABC);
    MyGraphType optiout;
    BranchDetection::SimplifyIt(ABC,optiout,15.0,"","pure");
    Write::GraphToVtk(f4,optiout);


//FlexibleComplex::AlphaTest(as,f3);
}

void UIHandler::HashComplexTest(std::list<Point> cloud)
{
    Alpha_shape_3 as(cloud.begin(),cloud.end());
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
    HashComplex s(cloud.size());
// Compute the time:
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    ComplexAlgorithms::LoadEmUp(as, s);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << "Loading up duration " << duration << std::endl;
    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
    int result = s.CrashTest();
    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
    std::cout << "Crashtest" << duration2 << std::endl;
}

void UIHandler::AlphaReebTest(std::list<Point> cloud, std::string outwrite, std::string inter)
{
    MyGraphType G;
    MyGraphType InterMediate;
    Computation::ComputeDeluanayTriangulation(G, cloud);
    double epsilon = 5.0;
    Computation::EpsilonSimplification(G, epsilon);
//Computation::BruteNeighborhoodGraph(G, cloud, epsilon);
    AlphaReeb_Parameters parameters(10.0,1.0);
    MyGraphType out;
    std::cout << "About to go into computation" << std::endl;
    AlphaReebComputation::Compute(G, parameters, out, InterMediate);
    Write::GraphToVtk(outwrite, out);
    Write::GraphToVtk(inter, InterMediate);



}

void UIHandler::MapperTest(std::string inwrite, std::string outwrite, Mapper_Parameters const & param)
{
    Filereader p;
    p.setPath(inwrite);
    std::list<Point> l;
    p.XYZRead(inwrite,l);
    std::vector<Point> v{ std::begin(l), std::end(l) };
    MyGraphType Final;
    Mapper ( v, param, Final);

    Write::GraphToVtk(outwrite, Final);
}

void UIHandler::ComputeAlphaReebInterior(MyGraphType & out, AlphaReeb_Parameters & param, std::vector<Point> & cloud, double epsilon)
{
    MyGraphType G;
    std::list<Point> cloudlist(cloud.begin(), cloud.end());
    Computation::ComputeDeluanayTriangulation(G, cloudlist);
    Computation::EpsilonSimplification(G, epsilon);
    MyGraphType Intermediate;
    AlphaReebComputation::Compute(G, param, out, Intermediate);

}


//void UIHandler::ComputeMapper(std::string inwrite, std::string outwrite)
//{
//    Filereader p;
//    p.setPath(inwrite);
//    std::list<Point> l;
//    p.XYZRead(inwrite,l);
//    std::vector<Point> v{ std::begin(l), std::end(l) };
//    MyGraphType Final;
//    Mapper ( v, param, Final);
//
//    Write::GraphToVtk(outwrite, Final);
//}

//        static void ComputeAlphaReeb(std::string inwrite, std::string outwrite, AlphaReeb_Parameters & param, double epsilon);


void UIHandler::ComputeAlphaReeb(std::string inwrite, std::string outwrite, AlphaReeb_Parameters & param, double epsilon)
{
    Filereader p;
    p.setPath(inwrite);
    std::list<Point> cloud;
    p.XYZRead(inwrite,cloud);
    MyGraphType G;
    MyGraphType InterMediate;
    Computation::ComputeDeluanayTriangulation(G, cloud);
    Computation::EpsilonSimplification(G, epsilon);
    MyGraphType out;
    AlphaReebComputation::Compute(G, param, out, InterMediate);
    Write::GraphToVtk(outwrite, out);
    //Write::GraphToVtk(inter, InterMediate);
}
void UIHandler::ThreeAlgorithmSuperTest(std::string inputfolder, std::string outputfolder, std::vector<std::string> inputs, std::vector<std::string> outputs)
{
    for (int i = 0; i < inputs.size(); i++)
    {
        for (int j = 0; j < 10; j++)
        {
            std::string numberj = std::to_string(j);
            std::string outputstring = outputfolder + numberj + "/";
            std::string inputstring = inputfolder + inputs[i];

            std::cout << "Starting round: " << numberj << ", outputting in: " << outputstring << std::endl;

            // Launch Standart:
            double bpar = 20 + 2*j;

            ClassicUIHandle smethod(bpar,1.5,1.5);
            double averageMST;
            std::string tmpstandart = outputstring + "Standart_" + outputs[i];
            smethod.CorrectCalculation(inputstring, tmpstandart, averageMST);

            std::cout << "Mapper time " << std::endl;
            // Launch Mapper:

            double parpar = 2 + 2*j;
            Mapper_Parameters param(15, 0.5, "Distance", parpar,parpar);
            tmpstandart = outputstring + "Mapper_" + outputs[i];
            UIHandler::MapperTest(inputstring, tmpstandart, param);

            // Launch AlphaReeb:

            std::cout << "AlphaReeeb time " << std::endl;
            double alpha = 2+2*j;
            AlphaReeb_Parameters parameters(alpha,1.0);
            tmpstandart = outputstring + "AlphaReeb_" + outputs[i];

            UIHandler::ComputeAlphaReeb(inputstring, tmpstandart, parameters, 3*averageMST);
        }
    }


}

// Parameters AlphaReeb: 10.0 20.0 30.0

void UIHandler::RunTestsForTheeAlgorithm(int n, int k, double minangle, double scale, int testN, double epsilon)
{
    //  std::vector<int> correctHomeo = {0,0,0};
//   std::vector<int> correctForm = {0,0,0};
    std::vector<std::string> filenames = {"normal.vtk", "mapper.vtk", "alphareeb.vtk"};
    for (int i = 0; i < testN; i++)
    {
        MyGraphType G;
        std::list<Point> cloud;
        GraphGeneration::RandomGraph1(k,minangle,scale,G);
        CloudGenerator::generatePoints(n, G, epsilon, cloud);
        ClassicUIHandle smethod(1.0,1.5,1.5);
        MyGraphType FG[3];
        double mstlength;
        MyGraphType mstree = Computation::computeMST(cloud);

        smethod.CorrectCalculationToGraph(FG[0], cloud, mstlength,"sd");

        Mapper_Parameters param(15, 0.5, "Distance", 5.0,5.0);
        AlphaReeb_Parameters AlphaParam(10.0, 1);
        double epsilon = 10*mstlength;
        std::vector<Point> vcloud{ std::begin(cloud), std::end(cloud)};
        Mapper ( vcloud, param, FG[1]);
        UIHandler::ComputeAlphaReebInterior(FG[2], AlphaParam, vcloud, epsilon);

        // Just print them:

        std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/output/Test5/" + std::to_string(i) + "/";
        std::cout << "Printing in folder: " << folder << std::endl;
        Write::GraphToVtk(folder + "mst.vtk", mstree);
        for (int pp = 0; pp < 3; pp++)
        {
            Write::GraphToVtk(folder + filenames[pp],FG[pp]);
            MyGraphType tmp;
            GraphSimplification::SimplifyGraph(FG[pp],tmp);
            bool correctForm = GraphSimplification::CorrectNumberOfEnds(k,tmp);
            bool correctHomeo = GraphSimplification::RecognizeStraGraph(tmp,k);
            double finalerror = Computation::AABBError(FG[pp],cloud);
            std::cout << "Round " << i << ": " << filenames[pp] << " | correct form: " << correctForm << " | correct homeo: " << correctHomeo << " | Error: " << finalerror << std::endl;


        }

        // Checking for corretness




    }
//    std::cout << "Number of correct recognitions: " << correct << std::endl;
}

void UIHandler::BetterTesting(int n, int k, double minangle, double scale, int testN, double epsilon)
{

    std::vector<std::string> filenames = {"normal.vtk", "mapper.vtk", "alphareeb.vtk"};
    MyGraphType G;
    std::list<Point> cloud;
    GraphGeneration::RandomGraph1(k,minangle,scale,G);
    CloudGenerator::generatePoints(n, G, epsilon, cloud);
    for (int i = 0; i < 10; i++)
    {

        ClassicUIHandle smethod(1.5,1.5,1.5);
        MyGraphType FG[3];
        double mstlength;
        MyGraphType mstree = Computation::computeMST(cloud);

        smethod.CorrectCalculationToGraph(FG[0], cloud, mstlength,"sd");

        Mapper_Parameters param(15, 0.5, "Distance", 0.5+1*i,0.5+1*i);
        AlphaReeb_Parameters AlphaParam(10.0 + i*5, 1);
        double epsilon = 10;
        std::vector<Point> vcloud{ std::begin(cloud), std::end(cloud)};
        Mapper ( vcloud, param, FG[1]);
        UIHandler::ComputeAlphaReebInterior(FG[2], AlphaParam, vcloud, epsilon);

        // Just print them:

        std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/output/Test9/" + std::to_string(i) + "/";
        std::cout << "Printing in folder: " << folder << std::endl;
        Write::GraphToVtk(folder + "mst.vtk", mstree);
        for (int pp = 0; pp < 3; pp++)
        {
            Write::GraphToVtk(folder + filenames[pp],FG[pp]);
            MyGraphType tmp;
            GraphSimplification::SimplifyGraph(FG[pp],tmp);
            bool correctForm = GraphSimplification::CorrectNumberOfEnds(k,tmp);
            bool correctHomeo = GraphSimplification::RecognizeStraGraph(tmp,k);
            double finalerror = Computation::AABBError(FG[pp],cloud);
            std::cout << "Round " << i << ": " << filenames[pp] << " | correct form: " << correctForm << " | correct homeo: " << correctHomeo << " | Error: " << finalerror << std::endl;


        }

        // Checking for corretness




    }
}

void UIHandler::FindBug(int n, int k, double minangle, double scale, int testN, double epsilon)
{
    MyGraphType output;
    MyGraphType G;
    std::list<Point> cloud;
    double mstlength;
    GraphGeneration::RandomGraph1(k,minangle,scale,G);
    CloudGenerator::generatePoints(n, G, epsilon, cloud);
    MyGraphType mstree = Computation::computeMST(cloud);
    Mapper_Parameters param(15, 0.5, "Distance", 3.0, 3.0);
    std::vector<Point> vcloud{ std::begin(cloud), std::end(cloud)};
    Mapper ( vcloud, param, output);
    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/output/Test8/";
    Write::GraphToVtk(folder + "mst.vtk", mstree);
    Write::GraphToVtk(folder + "output.vtk",output);



}


void UIHandler::DebugClassic(int n, int k, double minangle, double scale, int testN, double epsilon)
{
    for (int i = 0; i < testN; i++)
    {
        MyGraphType output;
        MyGraphType G;
        std::list<Point> cloud;
        double mstlength;
        GraphGeneration::RandomGraph1(k,minangle,scale,G);
        CloudGenerator::generatePoints(n, G, epsilon, cloud);
        MyGraphType mstree = Computation::computeMST(cloud);
        ClassicUIHandle smethod(50,1.5,4.0);
        smethod.CorrectCalculationToGraph(output, cloud, mstlength,"pure");
        std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/output/DebugTest/" + std::to_string(i) + "/";
        Write::GraphToVtk(folder + "mst.vtk", mstree);
        Write::GraphToVtk(folder + "output.vtk",output);


    }
}

void UIHandler::RunTheFinalTests(int k, double minangle, double scale, int testN, double epsilon, double alpha, double mappercluster)
{
    std::vector<int> successForm = {0,0,0};
    std::vector<int> successHomeo = {0,0,0};
    std::vector<double> avgTime = {0,0,0};
    std::vector<double> avgError = {0,0,0};
    std::vector<double> numberOfVertices = {0,0,0};

    for (int i = 0; i < testN; i++)
    {
        std::cout << "Round " << i << std::endl;
        double mstlength;
        MyGraphType FG[3];
        MyGraphType G;
        std::list<Point> cloud;
        GraphGeneration::RandomGraph1(k,minangle,scale,G);
        int n = 500*k;
        CloudGenerator::generatePoints(n, G, epsilon, cloud);
        std::vector<Point> vcloud{ std::begin(cloud), std::end(cloud)};

        ClassicUIHandle smethod(70.0,1.5,1.5);

        Mapper_Parameters param(15, 0.5, "Distance", mappercluster,mappercluster);
        AlphaReeb_Parameters AlphaParam(alpha, 1);

        for (int pp = 0; pp < 3; pp++)
        {
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            if (pp == 0)
            {
                smethod.CorrectCalculationToGraph(FG[0], cloud, mstlength,"pure");
            }
            else if (pp == 1)
            {
                Mapper ( vcloud, param, FG[1]);
            }
            else if (pp == 2)
            {
                UIHandler::ComputeAlphaReebInterior(FG[2], AlphaParam, vcloud, 10*mstlength);
            }
            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            MyGraphType tmp;
            GraphSimplification::SimplifyGraph(FG[pp],tmp);
            bool correctForm = GraphSimplification::CorrectNumberOfEnds(k,tmp);
            if (correctForm)
            {
                successForm[pp] = successForm[pp] + 1;
            }
            bool correctHomeo = GraphSimplification::RecognizeStraGraph(tmp,k);
            if (correctHomeo)
            {
                successHomeo[pp] = successHomeo[pp] + 1;
            }
            double finalerror = Computation::AABBError(FG[pp],cloud);
            std::cout << "pp: " << pp <<  " # Final error currently: " << finalerror << std::endl;
            avgError[pp] = avgError[pp] + finalerror;
            avgTime[pp] = avgTime[pp] + duration;
            numberOfVertices[pp] = numberOfVertices[pp] + boost::num_vertices(FG[pp]);
        }

    }



    std::vector<std::string> vectorStrings = {"Classic stats: ", "Mapper stats: ", "Alphareeb stats: "};
    for (int pp = 0; pp < 3; pp++)
    {
        double avgvalue = testN;
        avgTime[pp] = avgTime[pp] / avgvalue;
        avgError[pp] = avgError[pp] / avgvalue;
        numberOfVertices[pp] = numberOfVertices[pp] / avgvalue;
        std::cout << vectorStrings[pp] << std::endl;
        std::cout << "Success Form: " << successForm[pp] << std::endl;
        std::cout << "Success Homeo: " << successHomeo[pp] << std::endl;
        std::cout << "Average Time " << avgTime[pp] << std::endl;
        std::cout << "Average Error: " << avgError[pp] << std::endl;
        std::cout << "Average ammount of vertices: " << numberOfVertices[pp] << std::endl;
        std::cout << "---" << std::endl;
    }




}

void UIHandler::ComputeThingsRequiredForPaper()
{
    double mstlength = 0;
    std::string datafolder = "/home/yury/Dropbox/Projects/BranchPointFind/Data/";
    std::vector<std::string> FilePaths =  {"branched.xyz", "Christmas_present_end_tag.xyz", "worm.xyz"};
    std::vector<std::string> ThePaths(3);
    std::vector<std::string> fname = {"branched.vtk"}; //"christmas.vtk" , "worm.vtk"
    std::vector<std::string> ffname = {"branchedR.vtk"};
    std::vector<std::string> fffname = {"branchedRR.vtk"}; //,"christmasR.vtk" , "wormR.vtk"};
    for (int i = 0; i < ThePaths.size() ; i++)
    {
        ThePaths[i] = datafolder + FilePaths[i];
    }
    std::vector<double > parameters = {1.02};
    std::string outputfolder = "/home/yury/Dropbox/MicelleProject/Micelle/QuickPaperOutput/";

    for (int j = 0; j < parameters.size(); j++)
    {
        std::string numberj = std::to_string(j);
        std::string outf = outputfolder + numberj + "/";
        for (int i = 0; i < ffname.size(); i++)
        {
        MyGraphType outGraph;
       std::list<std::list<Point>> optipath;
       std::list<std::list<Point>> branchsimplified;
        Filereader p;
        p.setPath(ThePaths[i]);
        std::list<Point> cloud;
        p.XYZRead(ThePaths[i],cloud);
        ClassicUIHandle smethod(1.5,parameters[j],1.5);
        smethod.CorrectCalculationToGraphIntermediate(outGraph,optipath,branchsimplified, cloud, mstlength,"sd");
        Write::pathPrintToVtkPointlist(branchsimplified,outf + fffname[i]);
        Write::pathPrintToVtkPointlist(optipath,outf + ffname[i]);
        Write::GraphToVtk(outf + fname[i],outGraph);
        }

    }
}

//void UIHandler::PathSavior(std::list<Point> cloud, std::string f1, std::string f2, std::string f3, std::string f4, std::string f5)
//{
//Alpha_shape_3 as(cloud.begin(),cloud.end());
//Alpha_iterator opt = as.find_optimal_alpha(1);
//as.set_alpha(*opt);
//MyGraphType tree = Computation::computeMST(cloud);
//Write::GraphToVtk(f1,tree);
//Write::AlphaVTK(f2, as);
//MyGraphType ABC = FlexibleComplex::OptimizedSpanningTree(as, f3);
//Write::GraphToVtk(f3,ABC);
//MyGraphType optiout;
//BranchDetection::SimplifyIt(ABC,optiout,15.0,"");
//Write::GraphToVtk(f4,optiout);
//
//
////FlexibleComplex::AlphaTest(as,f3);
//}


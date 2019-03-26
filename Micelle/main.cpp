#include <iostream>
#include "Filereader.h"
#include <string>
#include <list>
#include <Definitions.h>
#include <Graph.h>
#include <Computation.h>
#include <Write.h>
#include "UIHandler.h"
#include <gudhi/Simplex_tree.h>
#include "AbstractComplex.h"
#include <boost/math/special_functions/binomial.hpp>
#include <chrono>
#include <queue>
#include <functional>
#include "HomologyComputator.h"
#include "binomial_coeff_table.h"
#include "HashForSimplex.h"
#include "ClassicUIHandle.h"
#include "DescriptorCluster.h"
#include "ClusterSort.h"
// ALPHA SHAPES:



void shoutOut(std::list<Point> points)
{
    for(auto it = points.begin(); it != points.end(); it++)
    {
        std::cout << (*it) << std::endl;
    }
}

void testminusone()
{

    std::string xyzfile = "/home/yury/Downloads/branched.xyz";
    Filereader p;
    p.setPath(xyzfile);
    std::list<Point> points;
    p.XYZRead(xyzfile,points);

    std::cout << "No bad alloc yet" << std::endl;

    UIHandler::HashComplexTest(points);

}

void testminusfive()
{
    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take23/";
    ClassicUIHandle b(50.0,1.5,1.5);
    b.GraphRecognitionTestSecond(4500,4,3,M_PI/3,100,10,10,folder);
}

void testminusfour()
{
    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take22/";
    ClassicUIHandle b(15000.0,1.5,1.5);
    b.GraphRecognitionTestFirst(4000,8,M_PI/3,100,10,20,folder);
}

void testminusthree()
{

    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take17/";
    ClassicUIHandle b(50.0,1.5,1.5);
    b.multipleTests(4000,8,M_PI/3,100,10,4, folder);

}

void testminustwo()
{
    std::cout << "ACCESSED " << std::endl;
    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take11/";
    std::string xyzfilez = "/home/yury/Downloads/branched.xyz";
    std::string f1 = folder + "Amst3.vtk";
    std::string f2 = folder + "Ashape3.vtk";
    std::string f3 = folder + "Optimst.vtk";
    std::string f4 = folder + "ApproxSkeleton.vtk";
    //std::string f5 = folder + "debug.txt";
    std::string f5 = folder + "StraightPath.vtk";
    ClassicUIHandle b(15.0,1.2,5.0);
    b.CalculateAll(xyzfilez,f1,f2,f3,f4,f5,"");
    std::cout << "Computation succeful" << std::endl;

}

void test()
{
    std::string folder = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take10/";
    std::string xyzfile = "/home/yury/Dropbox/Projects/BranchPointFind/Data/Christmas_present_end_tag.xyz";
    Filereader p;
    p.setPath(xyzfile);
    std::list<Point> points;
    // std::list<Point>* pp = & points;
    p.XYZRead(xyzfile,points);
    std::string f1 = folder + "Amst3.vtk";
    std::string f2 = folder + "Ashape3.vtk";
    std::string f3 = folder + "Optimst.vtk";
    std::string f4 = folder + "ApproxSkeleton.vtk";
    std::string f5 = folder + "debug.txt";
    UIHandler::AlphaDistanceComputer(points,f1,f2,f3,f4);
    std::cout << "Computation succeful" << std::endl;

}
void test2()
{
    std::cout << "Test beggining" << std::endl;
    UIHandler handler;
    handler.setFolder("/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/take12/");
    handler.multipleTests(4000,8,M_PI/3,100,10,4);
    std::cout << "Test succeful" << std::endl;

}

void printSimplex(ST & st, Simplex_handle & u)
{
    std::cout << "{ ";
    for (ST::Vertex_handle v : st.simplex_vertex_range(u))
    {
        std::cout << v << ", ";
    }
    std::cout << "}" << std::endl;
}

void test3()
{

    ST st;
    auto triangle012 = {0, 1, 2};
    auto triangle023 = {0, 2, 3};
    auto edge03 = {0, 3};
    st.insert_simplex_and_subfaces(triangle012);
    Simplex_handle kam = st.insert_simplex_and_subfaces(triangle023).first;
    st.insert_simplex_and_subfaces(edge03);
    auto edge02 = {0, 2};
    ST::Simplex_handle e = st.find(edge02);
    auto test = {3,0,2};
    ST::Simplex_handle testS = st.find(test);
    for (ST::Simplex_handle v : st.boundary_simplex_range(kam))
    {
        bool kumkum = v == e;
        std::string evaluation;
        if (kumkum)
        {
            evaluation = "true";
        }
        else
        {
            evaluation = "false";
        }
        printSimplex(st, v);
        std::cout << " Truth Value: " << evaluation << std::endl;
    }

    if (testS == kam)
    {
        std::cout << "The world is saved" << std::endl;
    }
    printSimplex(st,testS);


}

void test4()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    double abc = 10000;
    double k = 4;
    double result = boost::math::binomial_coefficient<double>(abc, k) +
                    boost::math::binomial_coefficient<double>(99999, 3) +
                    boost::math::binomial_coefficient<double>(99988, 2) +
                    9997;
    double result2 = result + result;

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << duration;
}

void test5()
{
    std::cout << "Starting test5" << std::endl;
    ST st;
    auto tetrah = {0,1,2,3,4};
    st.insert_simplex_and_subfaces(tetrah);
    //auto edge02 = {0, 2};
    std::vector<size_t> edge02 = {2,0};

    ST::Simplex_handle e = st.find(edge02);
    // ST::Simplex_handle e1337 = st.find(e);

    std::cout << "Dimension of the element:" << st.dimension(e) << std::endl;

    for (int i = 1; i < 5; i++)
    {
        std::cout << "Dimension " << i << std::endl;
        for (ST::Simplex_handle v : st.cofaces_simplex_range(e,i))
        {
            printSimplex(st,v);
        }

    }
    std::cout << "new one coming " << std::endl;
    auto edge03 = {0};
    ST::Simplex_handle sa = st.find(edge03);

    for (int i = 1; i < 5; i++)
    {
        std::cout << "Dimension " << i << std::endl;
        for (ST::Simplex_handle v : st.cofaces_simplex_range(sa,i))
        {
            printSimplex(st,v);
        }

    }


}
template<typename T> void print_queue(T& q)
{
    while(!q.empty())
    {
        std::cout << q.top() << " ";
        q.pop();
    }
    std::cout << '\n';
}
void test6()
{
    std::priority_queue<int> q;
    for(int n :
            {
                1,8,5,6,3,4,0,9,7,2
            })
    {
        q.push(n);
    }
    print_queue(q);
}

void test8()
{
    std::cout << "Starting test8" << std::endl;
    ST st;
    auto tetrah = {0,1,2,3,4};
    Simplex_handle kam = st.insert_simplex_and_subfaces(tetrah).first;
    for (auto v : st.boundary_simplex_range(kam))
    {
        printSimplex(st,v);
    }
// FIX THE PROBLEM USING HASHSET

}

struct Foo
{
    int k;
};

bool Compare(Foo a, Foo b)
{
    return a.k < b.k;
}

//typedef std::priority_queue<Foo, std::vector<Foo>, std::function<bool(Foo, Foo)>> customHeap;

void test7()
{
//customHeap pq(Compare);
    Foo one = {1};
    Foo two = {2};
    Foo three = {3};
    Foo four = {4};
    Foo five = {5};
    Foo six = {6};
    Foo seven  = {7};

    std::vector<Foo> w = {two,three,one,five,seven,six,four};
    std::make_heap(w.begin(), w.end(), Compare);
    for (auto it = w.begin(); it != w.end(); it++)
    {
        std::cout << (*it).k << std::endl;
    }

}

void test9()
{

    ST st;
    auto tetrah = {0,1,2,3,4};
    auto tetrah2 = {1,2,3,5,6};
    auto tetrah3 = {2,3,8};
    Simplex_handle kam = st.insert_simplex_and_subfaces(tetrah).first;
    Simplex_handle kam2 = st.insert_simplex_and_subfaces(tetrah2).first;
    st.insert_simplex_and_subfaces(tetrah3);
    auto simple = {1,2,3};
    Simplex_handle kam3 = st.find(simple);
    for(Simplex_handle v : st.star_simplex_range(kam3))
    {
        printSimplex(st,v);
    }


}

void test10()
{

    std::string xyzfile = "/home/yury/Downloads/worm.xyz";
    Filereader p;
    p.setPath(xyzfile);
    std::list<Point> points;
    // std::list<Point>* pp = & points;
    p.XYZRead(xyzfile,points);
    std::string f1 = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/Worm/Amst3.vtk";
    std::string f2 = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/Worm/Ashape3.vtk";
    UIHandler::AlphaSimplification(points, f1,f2);
    std::cout << "Computation succeful" << std::endl;



}

void sortTest()
{
    std::vector<int> w= {4,5,3,2,8,9,20,11};
    std::sort (w.rbegin(),w.rend());
    for (auto it = w.begin(); it != w.end(); it++)
    {
        std::cout << *it << " " ;
    }

}

void test20()
{
    Filereader p;
    std::string xyzfile = "/home/yury/Downloads/ring.xyz";
    std::list<Point> points;
    p.XYZRead(xyzfile,points);
    std::string f1 = "/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/Ring2/";
    Computation::FiltrationOfAlphaShapes(points, f1);



}

void test21()
{
    UIHandler handler;
    handler.setFolder("/home/yury/Dropbox/MicelleProject/Micelle/SuperOutput/");
    handler.multipleTestsPremium(10000,4,M_PI/8,100,10,10);



}

void BinomialHashTimeTest()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (int j = 0 ;  j < 100000; j++)
    {
        std::vector<size_t> t(4);
        int maxx = 10000;
        for (int i = 0; i < 4; i++)
        {
            int k = floor((maxx-(4-i))*Computation::unitRandom() + (4-i));
            t[i] = k;
            maxx = k;

        }
        Computation::BinomialHash(t);

    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << "Duration: " << duration << std::endl;


}

void testHomology()
{
    int coeff_field_characteristic = 1;
    Filtration_value min_persistence = 0.5;
    Simplex_tree st;
    auto triangle012 = {0, 1, 2};
    auto triangle013 = {0, 2, 3};
    auto triangle014 = {1, 2, 3};
    auto triangle015 = {0, 1, 3};
    st.insert_simplex_and_subfaces(triangle012, 1);
    st.insert_simplex_and_subfaces(triangle013, 1);
    st.insert_simplex_and_subfaces(triangle014, 1);
    st.insert_simplex_and_subfaces(triangle015, 1);

    Persistent_cohomology pcoh(st);
    pcoh.init_coefficients(coeff_field_characteristic);
//  pcoh.compute_persistent_cohomology(min_persistence);
// pcoh.output_diagram();
}

void testTriangle()
{
    Point a(-5,0,-5);
    Point b1(-10,1,-10);
    Point b2(-10,1,10);
    Point b3(10,1,10);
    Triangle k(b1,b2,b3);
//    std::cout << "Distance: " << sqrt(CGAL::squared_distance(k,a)) << std::endl;



}

void ActualHomology()
{
    HomologyComputator::TestRun();

}

void testBinomialCoefficent()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    binomial_coeff_table table(10000,4);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << "Duration: " << duration << std::endl;

}
void BinomialHashTimeTest2()
{
    HashForSimplex wasd = HashForSimplex(10000,4);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (int j = 0 ;  j < 100000; j++)
    {
        std::vector<int> t(4);
        int maxx = 10000;
        for (int i = 0; i < 4; i++)
        {
            int k = floor((maxx-(4-i))*Computation::unitRandom() + (4-i));
            t[i] = k;
            maxx = k;

        }
        wasd.giveMeHash(t);

    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << "Duration: " << duration << std::endl;
}

//typedef std::vector<int> wasdwasd;
//
//typedef std::unordered_map<wasdwasd, int, std::function<int(const wasdwasd&)>,
//std::function<bool(const wasdwasd&, const wasdwasd&)>> www;
//
//int hashFunction(const wasdwasd & a)
//{
//return a[0];
//}
//
//bool equals(wasdwasd & a, wasdwasd & b)
//{
//return true;
//}
//
//void hashsetTest()
//{
//www heap(hashFunction,equals);
//std::vector<int> abc = {2,1,0};
//heap.insert(std::make_pair(abc,1337));
//
//}

struct MyHash
{
    std::size_t operator()(const int& k) const
    {
        return 0;
    }
    bool operator()(const int& lhs, const int& rhs) const
    {
        return true;
    }


};

struct MyEqual
{
    bool operator()(const int& lhs, const int& rhs) const
    {
        return true;
    }
};

void hashTest()
{
    MyHash h;
    MyEqual e;
    std::unordered_map<int, std::string, MyHash, MyHash> m(42, h,h);
    m.insert(std::make_pair(1,"lol"));
    std::cout << m.find(1)->second << std::endl;

}

void unorderedSetTest()
{
    std::unordered_map<int,int> uorder;
    int bound = 100000;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < bound; i++)
    {
        uorder.insert(std::make_pair(i,i));
    }
    for (int i = bound; i < 2*bound; i++)
    {
        if (i % 10000 == 0)
        {
            std::cout << "Progresses" << std::endl;
            //    std::cout << uorder.find(i)->second << std::endl;
        }
    }


    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << "Duration: " << duration << std::endl;
}

void ReebGraphTest()
{
    std::string xyzfile = "/home/yury/Dropbox/Projects/BranchPointFind/Data/branched.xyz";
    Filereader p;
    p.setPath(xyzfile);
    std::list<Point> points;
    // std::list<Point>* pp = & points;
    p.XYZRead(xyzfile,points);
    std::string outwrite = "/home/yury/Dropbox/MicelleProject/Micelle/output/AlphaReeb/firstTest.vtk";
    std::string intermediate = "/home/yury/Dropbox/MicelleProject/Micelle/output/AlphaReeb/interTest.vtk ";
    UIHandler::AlphaReebTest(points, outwrite,intermediate);


}

void MapperTest()
{
    std::string xyzfile = "/home/yury/Dropbox/Projects/BranchPointFind/Data/branched.xyz";
    std::string outwrite = "/home/yury/Dropbox/MicelleProject/Micelle/output/Mapper/firstTest.vtk";
    Mapper_Parameters param(30, 0.25, "Distance", 1.0, 0.75);
    UIHandler::MapperTest(xyzfile, outwrite, param);
}

void AllocateCorrectly(std::vector<std::string> & model, std::vector<std::string> & input, std::vector<std::string> & output)
{
    for (int i = 0; i < model.size(); i++)
    {
        input[i] = model[i] + ".xyz";
        output[i] = model[i] + ".vtk";
    }
}


void AllThreeAlgorithmsOutput()
{
    std::string inputfolder = "/home/yury/Dropbox/Projects/BranchPointFind/Data/";
    std::vector<std::string> inputclouds = {"branched", "Christmas_present_end_tag", "worm", "multi_end_tag2"};
    std::string outputfolder = "/home/yury/Dropbox/MicelleProject/Micelle/output/AllThreeUpgrade/";
    std::vector<std::string> inputvector(inputclouds.size());
    std::vector<std::string> outputvector(inputclouds.size());
    AllocateCorrectly(inputclouds, inputvector, outputvector);
    UIHandler::ThreeAlgorithmSuperTest(inputfolder, outputfolder, inputvector, outputvector);
}

void debugclassicUi()
{
    double mstlength;
    std::string out = "/home/yury/Dropbox/MicelleProject/Micelle/output/debugBranched/theDebyg.vtk";
    std::string in = "/home/yury/Dropbox/Projects/BranchPointFind/Data/branched.xyz";
    ClassicUIHandle smethod(10,1.5,1.5);
    smethod.CorrectCalculation(in, out, mstlength);
    std::cout << "MSTLENGTHPARAMETER: " << mstlength << std::endl;
}
//4000,8,M_PI/3,100,10,4
void runTestsVeryNew()
{
    for (int i = 5; i < 6; i++)
    {
        UIHandler::RunTestsForTheeAlgorithm(4000, i, M_PI/3, 100, 1, 5);

    }
}

void anotheranotherTest()
{
    UIHandler::BetterTesting(1500, 5, M_PI/3, 100, 10, 3);
}

void DebugDebug()
{
    UIHandler::DebugClassic(4000, 5, M_PI/3, 100, 1, 5);

}

void sortTesting()
{
    std::vector<ClusterElement<vertex_descriptor>> descriptors;
    std::vector<vertex_descriptor> descriptorz= {vertex_descriptor(0),vertex_descriptor(1),2,3};
    std::vector<double> vv = {0.53, 0.12, 2.43, 1.43};
    for (int i = 0; i < 4; i++)
    {
        descriptors.push_back(ClusterElement<vertex_descriptor>(descriptorz[i],vv[i]));
    }


    std::sort(descriptors.rbegin(), descriptors.rend());

    for (auto it = descriptors.begin(); it != descriptors.end(); it++)
    {
        std::cout << (*it).returnObject() << " : " << (*it).returnValue() << " | ";
    }
    std::cout << std::endl;
    std::cout << "Succefully in this method " << std::endl;
}

void sortTestingTwo()
{
    std::vector<std::pair<vertex_descriptor,double>> descriptors;
    descriptors.push_back((std::make_pair(0,1)));
    descriptors.push_back((std::make_pair(1,0.67)));
    descriptors.push_back((std::make_pair(3,3.75)));
    descriptors.push_back((std::make_pair(2,4.5)));
//std::sort(descriptors.begin(), descriptors.end(),ClusterSort());

}

void DebugDebugDebug()
{

    UIHandler::FindBug(4000, 5, M_PI/3, 100, 1, 5);

}

void RunTheTests(int k, double alpha, double clustervalue)
{

    UIHandler::RunTheFinalTests(k,  M_PI/3, 100, 100, 5, alpha, clustervalue);


}

void actualProgramRun(int argc, char* argv[])
{

int kvalue = strtol(argv[1], NULL, 10);
double alpha =  atof(argv[2]);// 20;
double clustervalue = atof(argv[3]) ;//1.75;

std::cout << "Starting program with kvalue: " << kvalue << std::endl;
std::cout << "Starting program with alpha: " << alpha << std::endl;
std::cout << "Starting program with clustervalue: " << clustervalue << std::endl;

RunTheTests(kvalue, alpha, clustervalue);


}

int main(int argc, char* argv[])
{
    // std::cout << "WTF" << std::endl;
//   test2();
    // test5();
//  sortTest();
//test10();
//   test21();
// BinomialHashTimeTest();
//testTriangle();
//ActualHomology();
//test2();
//test2();
//testBinomialCoefficent();
//hashTest();
    //testminusone();
//BinomialHashTimeTest2();
//unorderedSetTest();
//test();
//testminustwo();
srand (time(NULL));
//ReebGraphTest();
//MapperTest();
//testminusfive();
//AllThreeAlgorithmsOutput();
//runTestsVeryNew();
//debugclassicUi();
//DebugDebug();
//sortTesting();
//sortTestingTwo();
//runTestsVeryNew();/
// anotheranotherTest();
// DebugDebugDebug();

//UIHandler::ComputeThingsRequiredForPaper();

//actualProgramRun(argc, argv);
anotheranotherTest();
std::cout << "Compiled succefully yes" << std::endl;


    return 0;
}




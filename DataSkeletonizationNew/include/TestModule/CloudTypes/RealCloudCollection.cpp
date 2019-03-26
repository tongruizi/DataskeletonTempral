#include "RealCloudCollection.h"
#include "AbstractCloudType.h"
#include "GeneralConvertor.h"
#include "SingleStar.h"
#include "DoubleStar.h"

RealCloudCollection::RealCloudCollection(std::string nameOfInstance, std::string pathName):
    AbstractCloudType(0,nameOfInstance)
{
    boost::filesystem::path p(pathName);
    std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(this->directories));
    this->number_of_runs = (this->directories).size();
    std::sort(this->directories.begin(), this->directories.end());
    this->correctGraphType.resize(this->number_of_runs);
    //(this->directories).size();
}
//    std::vector<*generatable> correctGraphType;

RealCloudCollection::~RealCloudCollection()
{
//! We destruct the pointers in the vector.
    for (auto it = this->correctGraphType.begin(); it != this->correctGraphType.end(); it++)
    {
        delete (*it);
    }
    this->correctGraphType.clear();
}


void RealCloudCollection::GenerateCloud(std::list<Point> & p, int iterationNumber)
{
    std::string kkk = (this->directories)[iterationNumber].path().string();
    GeneralConvertor::XYZtoPoint(kkk,p);
}


//static void RetriveGraphInformation(std::string filename, std::vector<int> & iterationanumber,
//std::vector<std::string> & graphType, std::vector<std::vector<int> parameters>);


//  virtual bool IsGraphCorrect(MyGraphType & G, int iterationnumber) =0;
// virtual bool DoesGraphHaveCorrectForm(MyGraphType & G, int iterationnumber) =0;

bool RealCloudCollection::IsGraphCorrect(MyGraphType & G, int iterationNumber)
{
    return this->correctGraphType[iterationNumber] -> IsGraphCorrect(G,iterationNumber);
}
bool RealCloudCollection::DoesGraphHaveCorrectForm(MyGraphType & G, int iterationNumber)
{
    return this->correctGraphType[iterationNumber] -> DoesGraphHaveCorrectForm(G,iterationNumber);

}

//! Unfourtunately for now, because C++ is not *something* languague, we have to use rude switch interface for this task.
void RealCloudCollection::SetCorrectnessOfGraphs(std::string filepath)
{
    std::vector<int> iterationanumber(this->number_of_runs, 0);
    std::vector<std::string> graphTypeString(this->number_of_runs,"");
    std::vector<std::vector<int>> parameters(this->number_of_runs,std::vector<int>());
    GeneralConvertor::RetriveGraphInformation(filepath,iterationanumber,graphTypeString,parameters);
    for (int i = 0; i < this->number_of_runs; i++)
    {
    std::cout << graphTypeString[i] << std::endl;
    }
    for (int i = 0; i < this->number_of_runs; i++)
    {
        std::string k = graphTypeString[i];
        generatable* g;

        if (k == "SingleStar")
        {
            g = new SingleStar(M_PI/3,parameters[i][0],1000,5,100,1,"SingleStar");
        }

        else if (k == "DoubleStar")
        {
            g = new DoubleStar(M_PI/3,parameters[i][0],parameters[i][1],1000,5,100,1,"DoubleStar");
        }
        else
        {
            std::cout << "Something went wrong in RealCloudCollection::SetCorrectnessOfGraphs method" << std::endl;
        }
        this->correctGraphType[i] = g;

    }
}
//GeneralConvertor::RetriveGraphInformation(filepath,)





//void FileSystemTest()
//{
////std::string p = "/home/yury/Dropbox/UnileverData/VTK_First_Files/";
//boost::filesystem::path p("/home/yury/Dropbox/UnileverData/VTK_First_Files/");
//std::vector<boost::filesystem::directory_entry> v; // To save the file names in a vector.
//
//    if(boost::filesystem::is_directory(p))
//    {
//        copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));
//        std::cout << p << " is a directory containing:\n";
//
//        for ( std::vector<boost::filesystem::directory_entry>::const_iterator it = v.begin(); it != v.end();  ++ it )
//        {
//            std::cout<< (*it).path().string()<<endl;
//        }
//    }
//
//}

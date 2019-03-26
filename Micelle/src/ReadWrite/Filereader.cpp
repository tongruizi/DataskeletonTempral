#include "Filereader.h"
#include "Definitions.h"
#include <fstream>
#include <iostream>

Filereader::Filereader()
{
    //ctor
}

void Filereader::setPath(std::string p)
{
    filepath = p;
}
std::string Filereader::returnPath()
{
    return filepath;
}

void Skipline(std::ifstream & s, int k)
{
    std::string es;
    for (int i = 0; i < k; i++)
    {
        getline(s,es);
    }

}

void Filereader::XYZRead(std::string ka, std::list<Point> & container)
{
    std::ifstream stream(ka.c_str());
    Point p;
    std::string s;
    Skipline(stream,2);
    while (stream >> s >> p)
    {
        container.push_back(p);
    }
}
void Filereader::customRead(std::string ka, std::list<Point> & container, int k)
{
    std::ifstream stream(ka.c_str());
    Point p;
    Skipline(stream,k);
    while (stream >> p)
    {
        container.push_back(p);
    }
}



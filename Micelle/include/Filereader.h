#ifndef FILEREADER_H
#define FILEREADER_H

#include <list>
#include <string>
#include <Definitions.h>


class Filereader
{
    std::string filepath;
    public:
        Filereader();
        void setPath(std::string path);
        std::string returnPath();
        void XYZRead(std::string ka, std::list<Point> & container);
        void customRead(std::string ka, std::list<Point> & container, int k);
    protected:

    private:
};

#endif // FILEREADER_H

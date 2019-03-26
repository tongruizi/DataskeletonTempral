#ifndef WRITE_H
#define WRITE_H
#include<Graph.h>
#include<Definitions.h>
#include<PointInfo.h>

class Write
{
    std::string folder;
    public:
        Write();
        static void GraphToVtk(std::string path, MyGraphType &G);
        static void AlphaVTK(const std::string &asfile, Alpha_shape_3 &as);
        static void BoundaryVerticesVTK(const std::string &asfile, Alpha_shape_3 &as);
        static void AlphaVTKSpecial(const std::string &asfile, std::vector<Chandler> & Simplices,
                            std::vector<Point> & verticesInfo, ST & st);
        static void SimplicesWriter(const std::string &asfile, std::vector<ST::Simplex_handle> handles[3], ST & st);
        static void DistancesFromBoundary(const std::string &asfile, std::map<Point, PointInfo> & dm);
        static void pathPrintToVtkPointlist(std::list<std::list<Point>> & paths, std::string directory);


    protected:

    private:
};

#endif // WRITE_H

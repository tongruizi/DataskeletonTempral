#include "GraphGeneration.h"
#include "Computation.h"

GraphGeneration::GraphGeneration()
{
    //ctor
}

double unitRandom()
{
    return ((double) rand() / (RAND_MAX));
}

double canAdd(std::list<Point> & points, Point q, Point origin, double minangle)
{

    double minAngle = 2*M_PI;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        Vector v1 = Vector(origin,q);
        Vector v2 = Vector(origin,(*it));
        // double angle =
        double cosine = v1 * v2 / CGAL::sqrt(v1*v1) / CGAL::sqrt(v2*v2);
        double angle = std::acos(cosine);
        // double angle = angleprime*2*M_PI/360;
        //    std::cout << angle << std::endl;
        minAngle = std::min(minAngle, angle);
    }
    return minAngle;
}

void GraphGeneration::RandomGraph1(int n, double minangle, double scale, MyGraphType & G)
{

    bool stillrun = true;
    double minangler;
    while(stillrun)
    {
        bool rerun = false;
        Point origin = Point(0,0,0);
        vertex_descriptor center = Graph::add_vertex(G,origin);
        std::list<Point> points;
        for (int i = 0; i < n; i++)
        {
            double t = unitRandom();
            double z = scale*t + (-scale)*(1-t);
            double u = unitRandom();
            double angl = u*2*M_PI;
            double y = scale*sqrt(1-(z*z)/(scale*scale))*sin(angl);
            double x = scale*sqrt(1-(z*z)/(scale*scale))*cos(angl);
            Point newP(x,y,z);
            minangler = canAdd(points, newP, origin,minangle);
            if(minangler < minangle)
            {
                rerun = true;
                break;
            }
            points.push_back(newP);
            vertex_descriptor boundary = Graph::add_vertex(G,newP);
            Graph::add_edge(G,center, boundary);
        }
        if (rerun)
        {
            G.clear();
            continue;
        }
        else
        {
            break;
        }
    }
}

vertex_descriptor GraphGeneration::TemplateFunction(Point extra, double xshift, double yshift, double zshift, MyGraphType & G,
int n, double minangle, double scale)
{
    bool stillrun = true;
    double minangler;
    Point origin = Point(0,0,0);
    std::list<Point> points;

    while(stillrun)
    {
        bool rerun = false;
        points.push_back(extra);
        for (int i = 0; i < n; i++)
        {
            double t = Computation::unitRandom();
            double z = scale*t + (-scale)*(1-t);
            double u = Computation::unitRandom();
            double angl = u*2*M_PI;
            double y = scale*sqrt(1-(z*z)/(scale*scale))*sin(angl);
            double x = scale*sqrt(1-(z*z)/(scale*scale))*cos(angl);
            Point newP(x,y,z);
            minangler = canAdd(points, newP, origin,minangle);
            if(minangler < minangle)
            {
                rerun = true;
                break;
            }
            points.push_back(newP);
          //  vertex_descriptor boundary = Graph::add_vertex(G,newP);
          //  Graph::add_edge(G,center, boundary);
        }
        if (rerun)
        {
            points.clear();
            continue;
        }
        else
        {
            break;
        }
    }
    Point neworigin(origin.x() + xshift, origin.y() + yshift, origin.z() + zshift);
    vertex_descriptor central = Graph::add_vertex(G,neworigin);
    for (auto it = points.begin(); it != points.end(); it++)
    {
    vertex_descriptor boundary = Graph::add_vertex(G,Point((*it).x() + xshift, (*it).y() + yshift, (*it).z() + zshift ));
    Graph::add_edge(G, central, boundary);
    }

    return central;

}

void GraphGeneration::RandomGraph2(int n, int n2, double minangle, double scale, MyGraphType & G)
{
vertex_descriptor first = GraphGeneration::TemplateFunction(Point(scale,0,0),0,0,0,G,n,minangle,scale);
vertex_descriptor second = GraphGeneration::TemplateFunction(Point(-scale,0,0),2*scale,0,0,G,n2,minangle,scale);
Graph::add_edge(G,first,second);

}



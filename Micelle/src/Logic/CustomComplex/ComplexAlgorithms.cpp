#include "ComplexAlgorithms.h"

//addVertex(int p)
ComplexAlgorithms::ComplexAlgorithms()
{
}

void ComplexAlgorithms::LoadEmUp(Alpha_shape_3 & as, HashComplex & h)
{

    std::map<Point, std::pair<size_t,bool>> points;
    size_t index = 0;
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        bool boundary = as.classify(v_it) == Alpha_shape_3::REGULAR;
        Point p = v_it -> point();
        points[p] = std::make_pair(index, boundary);
        index++;
    }
    std::cout << "Vertices loaded up " << std::endl;
    int j = 0;
    for (auto cell = as.finite_cells_begin(); cell != as.finite_cells_end(); cell++)
    {
    std::cout << "Handling cell" << j << std::endl;
        if (as.classify(cell) == Alpha_shape_3::INTERIOR)
        {
            Simplex q;
            bool boundary = false;
            for (int i = 0; i < 4; i++)
            {
                auto infopair = points.find(cell->vertex(0)->point())->second;
                q.push_back(infopair.first);
                if (infopair.second)
                {
                    boundary = true;
                }
            }
        std::sort(q.rbegin(),q.rend());
        h.addSimplex(q,boundary);
        }
        j++;
    }
}




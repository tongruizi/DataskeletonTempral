#include "Write.h"
#include "Graph.h"
#include "Definitions.h"


Write::Write()
{

}

void Write::GraphToVtk(std::string path, MyGraphType &G)
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
}

void Write::BoundaryVerticesVTK(const std::string &asfile, Alpha_shape_3 &as)
{
    int number = 0;
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        if (as.classify(v_it) != Alpha_shape_3::REGULAR)
        {
            number++;
        }
    }

    std::ofstream of(asfile);
    of << "# vtk DataFile Version 2.0\n\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";
    of << "POINTS " << number << " float\n";



    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        if (as.classify(v_it) != Alpha_shape_3::REGULAR)
        {
            of << v_it->point() << std::endl;
        }
    }


    of.close();


}

void Write::AlphaVTKSpecial(const std::string &asfile, std::vector<Chandler> & Simplices,
                            std::vector<Point> & verticesInfo, ST & st)
{
    int verticescount = verticesInfo.size();
    std::ofstream of(asfile);
    of << "# vtk DataFile Version 2.0\n\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";
    of << "POINTS " << verticescount << " float\n";
    int tetra_num = Simplices.size();
    int cellcount = tetra_num;
    for (int i = 0; i < verticescount; i++)
    {
        of << verticesInfo[i] << std::endl;
    }
    of << std::endl;
    of << "CELLS " << tetra_num << " " << 5 * tetra_num << std::endl;
    ST::Vertex_handle vH[4];
    for (auto it = Simplices.begin(); it != Simplices.end(); it++)
    {
        of << "4 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " " << (*it)[3] << std::endl;
    }
    of << std::endl;
    of << "CELL_TYPES " << tetra_num << std::endl;
    for (int i = 0; i < tetra_num; i++)
    {
        of << "10 ";
    }
    of << std::endl;
    of.close();



}

void Write::SimplicesWriter(const std::string &asfile, std::vector<ST::Simplex_handle> handles[3], ST & st)
{
    std::ofstream of(asfile);

    for (int i = 0; i < 3; i++)
    {
    std::cout << "Size of: " << i << handles[i].size() << std::endl;
       for (auto it = handles[i].begin(); it != handles[i].end(); it++)
        {
            std::cout << "Reached this statement" << std::endl;
            for (size_t u : st.simplex_vertex_range(*it))
            {
                of << u << " ";
                std::cout << "Vertex: " << u << std::endl;
            }
            of << std::endl;
        }


    }

    of.close();




}

void Write::DistancesFromBoundary(const std::string &asfile, std::map<Point, PointInfo> & dm)
{
    std::ofstream of(asfile);

    for (auto it = dm.begin(); it != dm.end(); it ++)
    {
    of << it -> first << " " << (it ->second).distance << std::endl;
    }

    of.close();


}

void Write::AlphaVTK(const std::string &asfile, Alpha_shape_3 &as)
{
    std::ofstream of(asfile);

// http://cgal-discuss.949826.n4.nabble.com/Help-with-filtration-and-filtration-with-alpha-values-td4659524.html#a4659549

    std::vector<Alpha_shape_3::Cell_handle> cells;
    std::vector<Alpha_shape_3::Facet> facets;
    std::vector<Alpha_shape_3::Edge> edges;
// tetrahedron = cell, they should be the interior, it is inside the 3D space
    as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::INTERIOR);
// triangles
// for the visualiization, don't need regular because tetrahedron will show it
//as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::SINGULAR);
// edges
    as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::SINGULAR);


    size_t tetra_num, tri_num, edge_num;
    tetra_num = cells.size();
    tri_num = facets.size();
    edge_num = edges.size();
// vertices: points <-> id
    std::map<Point, size_t> points;
    size_t index = 0;
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
// finite_.. is from DT class
        points[v_it->point()] = index;
        index++;
    }

// write
    of << "# vtk DataFile Version 2.0\n\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";
    of << "POINTS " << index << " float\n";
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        of << v_it->point() << std::endl;
    }

    of << std::endl;
    of << "CELLS " << tetra_num + tri_num + edge_num << " " << 5 * tetra_num + 4 * tri_num + 3 * edge_num << std::endl;
    for (auto cell:cells)
    {
        size_t v0 = points.find(cell->vertex(0)->point())->second;
        size_t v1 = points.find(cell->vertex(1)->point())->second;
        size_t v2 = points.find(cell->vertex(2)->point())->second;
        size_t v3 = points.find(cell->vertex(3)->point())->second;
        of << "4 " << v0 << " " << v1 << " " << v2 << " " << v3 << std::endl;
    }
// https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#ad6a20b45e66dfb690bfcdb8438e9fcae
    for (auto tri_it = facets.begin(); tri_it != facets.end(); ++tri_it)
    {
        of << "3 ";
        auto tmp_tetra = tri_it->first;
        for (int i = 0; i < 4; i++)
        {
            if (i != tri_it->second)
            {
                of << points.find(tmp_tetra->vertex(i)->point())->second << " ";
            }
        }
        of << std::endl;
    }
// https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#af31db7673a6d7d28c0bb90a3115ac695
    for (auto e : edges)
    {
        of << "2 ";
        auto tmp_tetra = e.get<0>();
        int p1, p2;
        p1 = e.get<1>();
        p2 = e.get<2>();
        of << points.find(tmp_tetra->vertex(p1)->point())->second << " "
           << points.find(tmp_tetra->vertex(p2)->point())->second << std::endl;
    }

    of << std::endl;
    of << "CELL_TYPES " << tetra_num + tri_num + edge_num << std::endl;
    for (int i = 0; i < tetra_num; i++)
    {
        of << "10 ";
    }
    for (int i = 0; i < tri_num; i++)
    {
        of << "5 ";
    }
    for (int i = 0; i < edge_num; i++)
    {
        of << "3 ";
    }
    of << std::endl;
    of.close();
}


void Write::pathPrintToVtkPointlist(std::list<std::list<Point>> & paths, std::string directory)
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
}







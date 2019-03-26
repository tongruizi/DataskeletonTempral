#include "AbstractComplex.h"
#include "Definitions.h"
#include "Computation.h"
#include "PointTools.h"
#include <queue>
#include <functional>
#include <utility>
#include <limits>
#include "Write.h"



AbstractComplex::AbstractComplex():
    sortguy(&(infoMap[3]))
{
    //  AbstractComplex::theComplex = ST();


}

void AbstractComplex::HandleToChandler(Simplex_handle & k, Chandler & p)
{
    for (size_t u : AbstractComplex::theComplex.simplex_vertex_range(k))
    {
        p.push_back(u);
    }
}

bool AbstractComplex::sizeOneAndCorrect(Simplex_handle & a, Simplex_handle & correct)
{
    int i = 1;
    bool detected = false;
    for (Simplex_handle v : AbstractComplex::theComplex.cofaces_simplex_range(a,1))
    {
        if ((i == 1) && (v == correct))
        {
            detected = true;
        }
        if (i == 2)
        {
            return false;
        }
        i++;
    }
    return detected;
}

bool AbstractComplex::freeBoundary(Simplex_handle & b)
{
    for (ST::Simplex_handle v : AbstractComplex::theComplex.boundary_simplex_range(b))
    {
        if (AbstractComplex::sizeOneAndCorrect(v, b))
        {
            return true;
        }
    }
    return false;
}

bool AbstractComplex::freeBoundarry(Chandler & www)
{
    Simplex_handle b = AbstractComplex::theComplex.find(www);
    return AbstractComplex::freeBoundary(b);
}

bool AbstractComplex::isFitC(int q,std::vector<ST::Vertex_handle> & b)
{
    Simplex_handle c = AbstractComplex::theComplex.find(b);
    return AbstractComplex::isFit(c,q);
}

bool AbstractComplex::isFit(Simplex_handle & b, int q)
{
    int op = 0;
    for(Simplex_handle v: AbstractComplex::theComplex.cofaces_simplex_range(b,1))
    {
        op++;
        if (op > q)
        {
            return false;
        }

    }
    return true;
}

int boolToInt(bool k)
{
    if (k == true)
    {
        return 1;
    }
    return 0;
}

void AbstractComplex::GrabVertices(Chandler & v, std::vector<Chandler> & h)
{
    for (int i = 0; i < 4; i++)
    {
        h.push_back({v[i]});
    }

}
void AbstractComplex::GrabEdges(Chandler & v, std::vector<Chandler> & h)
{
    h.push_back({v[0],v[1]});
    h.push_back({v[0],v[2]});
    h.push_back({v[0],v[3]});
    h.push_back({v[1],v[2]});
    h.push_back({v[1],v[3]});
    h.push_back({v[2],v[3]});

}

void AbstractComplex::GrabFaces(Chandler & v, std::vector<Chandler> & h)
{
    h.push_back({v[0],v[1],v[2]});
    h.push_back({v[0],v[1],v[3]});
    h.push_back({v[1],v[2],v[3]});
    h.push_back({v[0],v[2],v[3]});
}

int AbstractComplex::computeEuler(Chandler & v)
{
    int vertices = 0;
    int edges = 0;
    int faces = 0;
    ST::Vertex_handle v1[4][1];
    ST::Vertex_handle v2[6][2] = {{v[0],v[1]},{v[0],v[2]}, {v[0],v[3]},
        {v[1],v[2]}, {v[1],v[3]}, {v[2],v[3]}
    };
    ST::Vertex_handle v3[4][3] = {{v[0],v[1],v[2]}, {v[0],v[1],v[3]}, {v[1],v[2],v[3]}
        ,{v[0],v[2],v[3]}
    };
    v1[0][0] = v[0];
    v1[1][0] = v[1];
    v1[2][0] = v[2];
    v1[3][0] = v[3];
    for (size_t i = 0; i < 4; i++)
    {
        std::vector<ST::Vertex_handle> myvector(v1[i],v1[i] + + sizeof v1[i] /sizeof v1[i][0]);
        vertices = vertices + boolToInt(AbstractComplex::isFitC(3,myvector));
    }
    for (int i = 0; i < 6; i++)
    {
        std::vector<ST::Vertex_handle> myvector(v2[i],v2[i] + + sizeof v2[i] /sizeof v2[i][0]);
        edges = edges + boolToInt(AbstractComplex::isFitC(2,myvector));
    }
    for (int i = 0; i < 4; i++)
    {
        std::vector<ST::Vertex_handle> myvector(v3[i],v3[i] + + sizeof v3[i] /sizeof v3[i][0]);
        faces = faces + boolToInt(AbstractComplex::isFitC(1,myvector));
    }
    return vertices - edges + faces - 1;
}

bool AbstractComplex::AbleToDelete(Chandler & b)
{
    int k = AbstractComplex::computeEuler(b);
    return (k == 0);
    std::cout << "Euler number currently: " << k << std::endl;
}


//int AbstractComplex::findDistance(ST::Simplex_handle & a)
//{
//
//    ST* staddress = &(AbstractComplex::theComplex);
//    std::pair<ST*, ST::Simplex_handle> p = std::make_pair(staddress, a);
//    return (AbstractComplex::infoMap[3].find(p)->second).distance;
//
//
//}
//
//bool AbstractComplex::Compare(Simplex_handle & a, Simplex_handle & b)
//{
//    int length1 = AbstractComplex::findDistance(a);
//    int length2 = AbstractComplex::findDistance(b);
//    return length1 > length2;
//}

void AbstractComplex::addToInfoMap( Chandler & h,int dimension, bool boundary)
{
    std::sort(h.rbegin(), h.rend());
    CellInfo asdd = {boundary,false,false,std::numeric_limits<double>::max(),NULL};
    AbstractComplex::infoMap[dimension].insert(std::make_pair(h,asdd));
}

void AbstractComplex::writeSimplex(Simplex_handle & it, std::ofstream & of)
{
    for (size_t u : AbstractComplex::theComplex.simplex_vertex_range(it))
    {
        of << u << " ";
    }
    of << std::endl;

}

void AbstractComplex::Convert(Alpha_shape_3 & as)
{
    std::ofstream of("/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/Ring/debug2.txt");

    // std::vector<Alpha_shape_3::Cell_handle> cells;
    // std::vector<Alpha_shape_3::Facet> facets;
    //  std::vector<Alpha_shape_3::Edge> edges;
    // as.get_alpha_shape_cells(std::back_inserter(cells), Alpha_shape_3::INTERIOR);
    // as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::SINGULAR);
    // as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::SINGULAR);


    std::map<Point, size_t> points;
    size_t index = 0;
    of << "(bounddary) Vertices" << std::endl;
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++)
    {
        Point p = v_it -> point();
        points[p] = index;
        AbstractComplex::verticesInfo.push_back(p);
        Chandler apm = {index};
        ST::Simplex_handle h = AbstractComplex::theComplex.insert_simplex_and_subfaces(apm).first;
        bool boundary = as.classify(v_it) == Alpha_shape_3::REGULAR;
        AbstractComplex::addToInfoMap(apm,0,boundary);
        if(boundary)
        {
            AbstractComplex::boundaryContainer[0].push_back(apm);
            writeSimplex(h,of);
        }
        index++;
    }
    of << "(bounddary) Edges" << std::endl;


    for (auto edge = as.finite_edges_begin(); edge != as.finite_edges_end(); edge++)
    {
        bool k = as.classify(*edge) == Alpha_shape_3::EXTERIOR;
        if (k == false)
        {
            auto tmp_tetra = (*edge).get<0>();
            int p1, p2;
            p1 = (*edge).get<1>();
            p2 = (*edge).get<2>();
            size_t v0 = points.find(tmp_tetra->vertex(p1)->point())->second;
            size_t v1 = points.find(tmp_tetra->vertex(p2)->point())->second;
            Chandler apm = {v0,v1};
            bool boundary = as.classify(*edge) == Alpha_shape_3::REGULAR;
            ST::Simplex_handle h = AbstractComplex::theComplex.insert_simplex(apm).first;
            if(boundary)
            {
                AbstractComplex::boundaryContainer[1].push_back(apm);
                writeSimplex(h,of);

            }
            AbstractComplex::addToInfoMap(apm,1,boundary);
        }
    }


    of << "(bounddary) Faces" << std::endl;


    for (auto face = as.finite_facets_begin(); face != as.finite_facets_end(); face++)
    {
        bool k = as.classify(*face) == Alpha_shape_3::EXTERIOR;
        if (k == false)
        {
            auto tmp_tetra = face->first;
            Chandler apm(3);
            int j = 0;
            for (int i = 0; i < 4; i++)
            {
                if (i != face->second)
                {
                    apm[j] = points.find(tmp_tetra->vertex(i)->point())->second;
                    j++;
                }
            }
            bool boundary = as.classify(*face) == Alpha_shape_3::REGULAR;


            auto pair = AbstractComplex::theComplex.insert_simplex(apm);
            ST::Simplex_handle h = pair.first;
            bool istrue = pair.second ;
            if (istrue == true)
            {
                if(boundary)
                {
                    AbstractComplex::boundaryContainer[2].push_back(apm);
                    writeSimplex(h,of);

                }

                AbstractComplex::addToInfoMap(apm,2,boundary);
            }

        }
    }
    std::cout << "Facet Indexing completed" << std::endl;


    for (auto cell = as.finite_cells_begin(); cell != as.finite_cells_end(); cell++)
    {
        if (as.classify(cell) == Alpha_shape_3::INTERIOR)
        {
            size_t v0 = points.find(cell->vertex(0)->point())->second;
            size_t v1 = points.find(cell->vertex(1)->point())->second;
            size_t v2 = points.find(cell->vertex(2)->point())->second;
            size_t v3 = points.find(cell->vertex(3)->point())->second;
            Chandler apm = {v0,v1,v2,v3};
            auto paird = AbstractComplex::theComplex.insert_simplex(apm);
            if (paird.second == true)
            {
                AbstractComplex::addToInfoMap(apm,3,false);
            }
        }
    }
    std::cout << "Cell Indexing completed" << std::endl;

    of.close();
}


typedef std::priority_queue<Simplex_handle, std::vector<Simplex_handle>,
        std::function<bool(Simplex_handle, Simplex_handle)>> simplexHeap;


Point AbstractComplex::BaryCenterOfSimplex(Simplex_handle h)
{
    std::vector<Point> p;
    for (size_t qp : AbstractComplex::theComplex.simplex_vertex_range(h))
    {
        Point q = AbstractComplex::verticesInfo[qp];
        p.push_back(q);
    }
    return PointTools::findBarycenter(p);
}

Point AbstractComplex::BaryCenterOfSimplex(Chandler & v)
{
    std::vector<Point> p;
    for (int i = 0; i < v.size(); i++)
    {
        Point q = AbstractComplex::verticesInfo[v[i]];
        p.push_back(q);
    }
    return PointTools::findBarycenter(p);
}

void AbstractComplex::MarkBoundary()
{
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < AbstractComplex::boundaryContainer[j].size(); i++)
        {
            std::cout << "J: " << j << "I: " << i << std::endl;
            Chandler e = boundaryContainer[j][i];
            Point a = AbstractComplex::BaryCenterOfSimplex(e);
            Simplex_handle h = AbstractComplex::theComplex.find(e);
            for (Simplex_handle qq : AbstractComplex::theComplex.cofaces_simplex_range(h, 3-j))
            {
                Chandler nrr;
                AbstractComplex::HandleToChandler(qq, nrr);
                Point b = AbstractComplex::BaryCenterOfSimplex(nrr);
                double newd = sqrt(CGAL::squared_distance(a,b));
                double oldd =  (AbstractComplex::infoMap[3].find(nrr)->second).distance;
                if (newd < oldd)
                {
                    (AbstractComplex::infoMap[3].find(nrr)->second).distance = newd;
                    (AbstractComplex::infoMap[3].find(nrr)->second).prev = &e;
                }
            }
        }

    }
}

bool AbstractComplex::handleSimplex(Chandler & to, Chandler & from,
                                    Chandler &mid, std::vector<Chandler> & heap)

{
//  Chandler hp;
    // AbstractComplex::HandleToChandler(to, hp);
    //  Chandler hp2;
//   AbstractComplex::HandleToChandler(from, hp2);


    bool istrue = false;
    Point topoint = AbstractComplex::BaryCenterOfSimplex(to);
    Point frompoint = AbstractComplex::BaryCenterOfSimplex(from);
    Point midpoint = AbstractComplex::BaryCenterOfSimplex(mid);
    double distanceToBase = (AbstractComplex::infoMap[3].find(from)->second).distance;
    double newdistance = sqrt(CGAL::squared_distance(topoint, midpoint)) +
                         sqrt(CGAL::squared_distance(frompoint, midpoint)) + distanceToBase;
    double olddistance = (AbstractComplex::infoMap[3].find(to)->second).distance;

    // Find the corresponding sequence element:


    if ((AbstractComplex::infoMap[3].find(to)->second).boundary == true)
    {
        if (newdistance < olddistance)
        {
            (AbstractComplex::infoMap[3].find(to)->second).distance = newdistance;
            (AbstractComplex::infoMap[3].find(to)->second).prev = NULL;
            istrue = true;
        }
    }
    else
    {
        if (newdistance < olddistance)
        {
            (AbstractComplex::infoMap[3].find(to)->second).distance = newdistance;
            (AbstractComplex::infoMap[3].find(to)->second).prev = NULL;
        }
        if(AbstractComplex::freeBoundarry(from))
        {
            (AbstractComplex::infoMap[3].find(to)->second).boundary = true;
            heap.push_back(to);
            std::push_heap(heap.begin(),heap.end(),AbstractComplex::sortguy);
        }

    }

    return istrue;
}


bool AbstractComplex::Checker(std::vector<Chandler> & faces, int p, Simplex_handle & topHandle,
                              std::vector<Chandler> & heap, std::vector<Chandler> & deletionList, Chandler & top
                             )
{

//  Chandler btc;
    // AbstractComplex::HandleToChandler(to, btc);
    //  Chandler gtc;
//   AbstractComplex::HandleToChandler(from, gtc);


    bool needToRebuild = false;
    for (auto it = faces.begin(); it != faces.end(); it++)
    {
        ST::Simplex_handle gt = AbstractComplex::theComplex.find(*it);
        int cofaces = 0;
        for (ST::Simplex_handle bt : AbstractComplex::theComplex.cofaces_simplex_range(gt,p))
        {
            cofaces++;
            if (bt != topHandle)
            {
                // bool AbstractComplex::handleSimplex(Simplex_handle & to, Simplex_handle & from,
                //  Simplex_handle &mid, std::vector<Simplex_handle> & heap)
                Chandler btc;
                AbstractComplex::HandleToChandler(bt, btc);
                Chandler gtc;
                AbstractComplex::HandleToChandler(gt, gtc);
                bool k = AbstractComplex::handleSimplex(btc, top,gtc, heap);

                if (k)
                {
                    needToRebuild = true;
                }
            }
        }
        if(cofaces == 0)
        {
            deletionList.push_back(*it);
        }
    }

    return needToRebuild;


}

void AbstractComplex::SimplifyComplex()
{

// DEBUG:
    std::cout << "Debugging" << std::endl;
//   Write::SimplicesWriter("/home/yury/Dropbox/MicelleProject/Micelle/ActualOutput/Ring/debug.txt", AbstractComplex::boundaryContainer, AbstractComplex::theComplex);


//

    std::vector<Chandler> heap;
    ST* address = &(AbstractComplex::theComplex);
// MarkBoundary;
    std::cout << "Before mark boundary" << std::endl;
    AbstractComplex::MarkBoundary();
    std::cout << "After mark boundary" << std::endl;
    int jjj =0;
    for (auto it = AbstractComplex::infoMap[3].begin(); it != AbstractComplex::infoMap[3].end(); it++)
    {
        Simplex_handle b = AbstractComplex::theComplex.find((*it).first);
        std::cout << "Round " << jjj << std::endl;
        if (AbstractComplex::freeBoundary(b))
        {
            heap.push_back((*it).first);
            (*it).second.boundary = true;
        }
        jjj++;
    }
    std::make_heap(heap.begin(), heap.end(), AbstractComplex::sortguy);
    std::cout << "LOL" << std::endl;
    std::cout << "HEAPSIZE CURRENTLY: " << heap.size() << std::endl;
    while(heap.size() > 0)
    {
        std::vector<Chandler> deletionList;
        std::cout << "Heapsize: " << heap.size() << std::endl;
        Chandler top = heap.front();
        deletionList.push_back(top);
        Simplex_handle topHandle = AbstractComplex::theComplex.find(top);
        std::pop_heap(heap.begin(), heap.end(), AbstractComplex::sortguy);
        heap.pop_back();
        if (AbstractComplex::AbleToDelete(top))
        {
            bool needToRebuild = false;
            std::vector<Chandler> faces;
            std::vector<Chandler> edges;
            std::vector<Chandler> verticess;
            AbstractComplex::GrabFaces(top, faces);
            AbstractComplex::GrabEdges(top, edges);
            AbstractComplex::GrabVertices(top, verticess);

            // Roll Over Facets:

            //bool AbstractComplex::Checker(std::vector<Chandler> & faces,
            // int p, Simplex_handle & topHandle,
            //std::vector<Simplex_handle> & heap

            bool k0 = AbstractComplex::Checker(faces, 1, topHandle, heap,deletionList, top);
            // Roll Over Edges:

            bool k1 = AbstractComplex::Checker(edges, 2, topHandle, heap,deletionList, top);
            // Roll Over Vertices:
            bool k2 = AbstractComplex::Checker(verticess, 3, topHandle, heap,deletionList, top);

            if (k0 || k1 || k2)
            {
                needToRebuild = true;
            }

            if (needToRebuild)
            {
                std::make_heap(heap.begin(), heap.end(),AbstractComplex::sortguy);
            }

            for (int i = 0; i < deletionList.size(); i++)
            {
                AbstractComplex::theComplex.remove_maximal_simplex(AbstractComplex::theComplex.find(deletionList[i]));
            }
        }
        else
        {
            //output.insert_simplex_and_subfaces(top);
            (AbstractComplex::infoMap[3].find(top)->second).irreducable = true;
        }

    }
}

void AbstractComplex::CollectIrreducableCells(std::vector<Chandler> & simplices)
{
    for (auto it = AbstractComplex::infoMap[3].begin(); it != AbstractComplex::infoMap[3].end(); it++)
    {
        if ((*it).second.irreducable == true)
        {
            simplices.push_back((*it).first);
        }

    }


}

ST* AbstractComplex::returnST()
{
    return &(AbstractComplex::theComplex);

}
std::vector<Point>* AbstractComplex::returnVertexIndexation()
{
    return &(AbstractComplex::verticesInfo);

}


// Boundary marked:






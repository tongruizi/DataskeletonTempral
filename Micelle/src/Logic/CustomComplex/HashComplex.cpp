#include "HashComplex.h"
#include <iostream>
#include <chrono>
#include <ctime>



// typedef std::unordered_set<Simplex, HashForSimplex, EqualsComplex> SimplexSet;


HashComplex::HashComplex(int n):
    simplexHash(n,4), boundaryCells(10,simplexHash, simplexEquals), interiorCells(10,simplexHash, simplexEquals)
{


}

void interSection(Simplex & intersection, Simplex & k1, Simplex & k2)
{
    int j = 0;
    int i = 0;

    while ((i < k1.size()) && (j < k2.size()))
    {
        if (k1[i] == k2[j])
        {
            intersection.push_back(k1[i]);
        }
    }
    if ((k1[i] > k2[j]))
    {
        i++;
    }
    else
    {
        j++;
    }
}

int HashComplex::EulerChangeComputator(Simplex & s, SimplexCount & m)
{
    SimplexSet w[3] = {SimplexSet(10,simplexHash,simplexEquals),SimplexSet(10,simplexHash,simplexEquals),
    SimplexSet(10,simplexHash,simplexEquals)};

    for (auto it = m.begin(); it != m.end(); it++)
    {
        Simplex intersection;
        Simplex k1 = it -> first;
        interSection(intersection, k1, s);
        w[intersection.size()-1].insert(intersection);
    }

    return -1+(4-w[2].size()) - (6-w[1].size()) + (4-w[0].size());
}

void HashComplex::getNeighbors(Simplex & s, SimplexCount & m)
{
    for (int i = 0; i < s.size(); i++)
    {
        for (auto it = (vertices.find(s[i])->second).cells.begin();
                it != (vertices.find(s[i])->second).cells.end(); it++)
        {
            auto git = m.find(*it);
            if (git == m.end())
            {
                m.insert(std::make_pair(*it,1));
            }
            else
            {
                (git -> second) = (git -> second) + 1;
            }
        }
    }
}

void HashComplex::removeSimplex(Simplex & s, bool boundary)
{
    for (int i = 0; i < s.size(); i++)
    {
        auto found = vertices.find(s[i]);
        (found->second).cells.erase(s);
    }

    if (boundary)
    {
        boundaryCells.erase(s);
    }
    else
    {
        interiorCells.erase(s);
    }
}

void HashComplex::addSimplex(Simplex & s, bool boundary)
{
std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < s.size(); i++)
    {
        (vertices.find(s[i])->second).cells.insert(s);
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << "Duration Simplex additions to vertices: " << duration << std::endl;

    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

    if (boundary)
    {
        boundaryCells.insert(s);
    }
    else
    {
        interiorCells.insert(s);
    }
std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();

    std::cout << "Duration Simplex additions step2: " << duration2 << std::endl;
}

int HashComplex::CrashTest()
{
    int totalcount = 0;
    for (auto it = boundaryCells.begin(); it != boundaryCells.end(); it++)
    {
        SimplexCount m(10,simplexHash,simplexEquals);
        auto deref = *it;
        HashComplex::getNeighbors(deref, m);
        int op = HashComplex::EulerChangeComputator(deref, m);
        if (op == 0)
        {
            totalcount++;
        }
    }
    return totalcount;
}

void HashComplex::addVertex(int p)
{
vertices.insert(std::make_pair(p,Vertex{p,SimplexSet(42,simplexHash, simplexEquals)}));
}

//    for (int i = 0; i < 4; i++)
//    {
//        w[0].insert(std::vector<int> {s[i]});
//    }
//    for (int i = 0; i < 3; i++)
//    {
//        for (int j = i+1; j < 4)
//        {
//            w[1].insert(std::vector<int> {s[i],s[j]});
//        }
//    }
//    for (int i = 0; i < 4; i++)
//    {
//        std::vector<int> p;
//        for (int j = 0; j < 4; j++)
//        {
//            if (i != j)
//            {
//                p.push_back(s[j]);
//            }
//        }
//        w[2].insert(p);
//    }

// loop over simpölices:


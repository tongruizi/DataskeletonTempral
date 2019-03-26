#ifndef HASHCOMPLEX_H
#define HASHCOMPLEX_H
#include "Definitions.h"
#include "HashForSimplex.h"
#include "EqualsComplex.h"
struct Vertex;
struct Cell;

typedef std::vector<int> Simplex;
typedef std::unordered_map<Simplex, int, HashForSimplex, EqualsComplex> SimplexCount;
typedef std::unordered_set<Simplex, HashForSimplex, EqualsComplex> SimplexSet;
typedef std::unordered_map<int,Vertex> vertexInfo;

//
struct Vertex
{
int id;
SimplexSet cells;
};
//struct Cell
//{
//double distance;
//Point barycenter;
//};


class HashComplex
{
HashForSimplex simplexHash;
SimplexSet boundaryCells;
SimplexSet interiorCells;
//SimplexHash Cells;
vertexInfo vertices;
EqualsComplex simplexEquals;

    public:
        HashComplex(int a);
        void addSimplex(Simplex & s, bool boundary);
        int CrashTest();
        void removeSimplex(Simplex & s, bool boundary);
        void getNeighbors(Simplex & s, SimplexCount & m);
        int EulerChangeComputator(Simplex & s, SimplexCount & m);
        void addVertex(int p);


    protected:

    private:
};

#endif // HASHCOMPLEX_H

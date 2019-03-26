#ifndef ABSTRACTCOMPLEX_H
#define ABSTRACTCOMPLEX_H

#include <Definitions.h>
#include <unordered_map>
#include "SimplexHash.h"
#include "SimplexComparator.h"
#include "SimplexSort.h"
#include "CellInfo.h"

class AbstractComplex
{
ST theComplex;
std::vector<Chandler> boundaryContainer[3];
std::vector<Point> verticesInfo;
std::unordered_map<Chandler, CellInfo,
SimplexHash, SimplexComparator> infoMap[4];
SimplexSort sortguy;
public:
    AbstractComplex();
    void HandleToChandler(Simplex_handle & k, Chandler & p);
    void Convert(Alpha_shape_3 & as);
    void addToInfoMap( Chandler & h,int dimension, bool boundary);
   // bool Compare(Simplex_handle & a, Simplex_handle & b);
  //  int  findDistance(Simplex_handle & a);
    Point BaryCenterOfSimplex(Simplex_handle h);
    Point BaryCenterOfSimplex(Chandler & v);
    void SimplifyComplex();
    bool AbleToDelete(ST::Simplex_handle & b);
    bool AbleToDelete(Chandler & b);
    void MarkBoundary();
    bool sizeOneAndCorrect(Simplex_handle & a, Simplex_handle & correct);
    bool freeBoundarry(Chandler & www);
    bool freeBoundary(Simplex_handle & b);
    bool isFitC(int q, std::vector<ST::Vertex_handle> & b);
    bool isFit(Simplex_handle & b, int q);
    int computeEuler(Chandler & b);
    void GrabVertices(Chandler & v, std::vector<Chandler> & h);
    void GrabEdges(Chandler & v, std::vector<Chandler> & h);
    void GrabFaces(Chandler & v, std::vector<Chandler> & h);
    bool handleSimplex(Chandler & to, Chandler & from,
                                    Chandler &mid, std::vector<Chandler> & heap);
    void CollectIrreducableCells(std::vector<Chandler> & simplices);
    ST* returnST();
    std::vector<Point>* returnVertexIndexation();
    bool Checker(std::vector<Chandler> & faces, int p, Simplex_handle & topHandle,
    std::vector<Chandler> & heap, std::vector<Chandler> & deletionList, Chandler & top);

protected:

private:
void writeSimplex(Simplex_handle & it, std::ofstream & of);

};

#endif // ABSTRACTCOMPLEX_H

#ifndef HOMOLOGYCOMPUTATOR_H
#define HOMOLOGYCOMPUTATOR_H


#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <iostream>
#include <ctime>
#include <utility>
#include <vector>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;
using typeVectorVertex = std::vector< Simplex_tree::Vertex_handle >;

class HomologyComputator
{
    public:
        HomologyComputator();
        static void TestRun();
        static void Compute(Simplex_tree & st);

    protected:

    private:
};

#endif // HOMOLOGYCOMPUTATOR_H

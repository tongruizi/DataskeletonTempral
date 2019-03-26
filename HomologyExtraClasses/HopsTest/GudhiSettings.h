#ifndef GUDHISETTINGS_H
#define GUDHISETTINGS_H

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <iostream>
#include <ctime>
#include <utility>
#include <vector>

//struct MyOptions : Gudhi::Simplex_tree_options_full_featured {
//  // Not doing persistence, so we don't need those
//  static const bool store_key = false;
//  static const bool store_filtration = false;
//  // I have few vertices
//  typedef short Vertex_handle;
//};

typedef Gudhi::Simplex_tree<> ST;
typedef ST::Filtration_value Filtration_value;
typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
typedef Gudhi::persistent_cohomology::Persistent_cohomology<ST, Field_Zp > Persistent_cohomology;
typedef std::vector< ST::Vertex_handle > typeVectorVertex;

#endif // GUDHISETTINGS_H

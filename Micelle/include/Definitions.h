#ifndef DEFINITONS_H
#define DEFINITONS_H

#include <queue>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/program_options.hpp>
#include <boost/graph/graphviz.hpp>
#include <utility>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <unordered_map>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/graph_simplicial_complex.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>



typedef std::vector<size_t> Chandler;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> Triangulation;
typedef Kernel::Point_3 Point;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Triangle_3 Triangle;

typedef CGAL::Alpha_shape_vertex_base_3<Kernel>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>      Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>             Alpha_shape_3;
typedef Alpha_shape_3::Alpha_iterator                    Alpha_iterator;

typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef std::list<Segment>::iterator IteratorSegment;
typedef CGAL::AABB_segment_primitive<Kernel, IteratorSegment> PrimitiveSegment;
typedef CGAL::AABB_traits<Kernel, PrimitiveSegment> TraitsSegment;
typedef CGAL::AABB_tree<TraitsSegment> SegmentTree;

struct MyOptions : Gudhi::Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static const bool store_key = false;
  static const bool store_filtration = false;
  // I have few vertices
  typedef size_t vertex_handle;
};

typedef Gudhi::Simplex_tree<MyOptions> ST;
typedef ST::Simplex_handle Simplex_handle;
typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef Simplex_tree::Filtration_value Filtration_value ;
typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp ;
typedef Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp > Persistent_cohomology;
#endif



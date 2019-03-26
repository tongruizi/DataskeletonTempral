/**
 * @file traits.hpp
 * @author Ryan Curtin
 *
 * Specialization of the TreeTraits class for the BinarySpaceTree type of tree.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef MLPACK_CORE_TREE_BINARY_SPACE_TREE_TRAITS_HPP
#define MLPACK_CORE_TREE_BINARY_SPACE_TREE_TRAITS_HPP

#include <mlpack/core/tree/tree_traits.hpp>
#include <mlpack/core/tree/ballbound.hpp>

namespace mlpack {
namespace tree {

/**
 * This is a specialization of the TreeTraits class to the BinarySpaceTree tree
 * type.  It defines characteristics of the binary space tree, and is used to
 * help write tree-independent (but still optimized) tree-based algorithms.  See
 * mlpack/core/tree/tree_traits.hpp for more information.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
             class SplitType>
class TreeTraits<BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
                                 SplitType>>
{
 public:
  /**
   * Each binary space tree node has two children which represent
   * non-overlapping subsets of the space which the node represents.  Therefore,
   * children are not overlapping.
   */
  static const bool HasOverlappingChildren = false;

  /**
   * Some segments are split into multiple cells
   */
  static const bool HasDuplicatedPoints = true;

  /**
   * There is no guarantee that the first point in a node is its centroid.
   */
  static const bool FirstPointIsCentroid = false;

  /**
   * Points are not contained at multiple levels of the binary space tree.
   */
  static const bool HasSelfChildren = false;

  /**
   * Points are not rearranged during building of the tree.
   */
  static const bool RearrangesDataset = false;

  /**
   * This is always a binary tree.
   */
  static const bool BinaryTree = true;

  /**
   * There are dupliacted points, so the number of descendants is not unique
   *
   */
  static const bool UniqueNumDescendants = false;
};


} // namespace tree
} // namespace mlpack

#endif

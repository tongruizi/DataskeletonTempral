/**
 * @file binary_space_tree.hpp
 *
 * Definition of generalized binary space partitioning tree (BinarySpaceTree).
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef MLPACK_CORE_TREE_BINARY_SPACE_SEGMENT_TREE_BINARY_SPACE_TREE_HPP
#define MLPACK_CORE_TREE_BINARY_SPACE_SEGMENT_TREE_BINARY_SPACE_TREE_HPP

#include <mlpack/prereqs.hpp>

#include "../statistic.hpp"
#include "midpoint_split.hpp"

namespace mlpack {
namespace tree /** Trees and tree-building procedures. */ {

/**
 * A binary space partitioning tree, such as a KD-tree or a ball tree.  Once the
 * bound and type of dataset is defined, the tree will construct itself.  Call
 * the constructor with the dataset to build the tree on, and the entire tree
 * will be built.
 *
 * This particular tree does not allow growth, so you cannot add or delete nodes
 * from it.  If you need to add or delete a node, the better procedure is to
 * rebuild the tree entirely.
 *
 * This tree does take one runtime parameter in the constructor, which is the
 * max leaf size to be used.
 *
 * @tparam MetricType The metric used for tree-building.  The BoundType may
 *     place restrictions on the metrics that can be used.
 * @tparam StatisticType Extra data contained in the node.  See statistic.hpp
 *     for the necessary skeleton interface.
 * @tparam MatType The dataset class.
 * @tparam BoundType The bound used for each node.  HRectBound, the default,
 *     requires that an LMetric<> is used for MetricType (so, EuclideanDistance,
 *     ManhattanDistance, etc.).
 * @tparam SplitType The class that partitions the dataset/points at a
 *     particular node into two parts. Its definition decides the way this split
 *     is done.
 */
template<typename MetricType,
         typename StatisticType = EmptyStatistic,
         typename MatType = arma::mat,
         template<typename BoundMetricType, typename...> class BoundType =
            bound::HRectBound,
         template<typename SplitBoundType, typename SplitMatType>
            class SplitType = MidpointSplit>
class BinarySpaceSegmentTree
{
 public:
  //! So other classes can use TreeType::Mat.
  typedef MatType Mat;
  //! The type of element held in MatType.
  typedef typename MatType::elem_type ElemType;

  typedef SplitType<BoundType<MetricType>, MatType> Split;

 private:
  //! The left child node.
  BinarySpaceSegmentTree* left;
  //! The right child node.
  BinarySpaceSegmentTree* right;
  //! The parent node (NULL if this is the root of the tree).
  BinarySpaceSegmentTree* parent;
  //! The bound object for this node.
  BoundType<MetricType> bound;
  //! Any extra data contained in the node.
  StatisticType stat;
  //! The distance from the centroid of this node to the centroid of the parent.
  ElemType parentDistance;
  //! The worst possible distance to the furthest descendant, cached to speed
  //! things up.
  ElemType furthestDescendantDistance;
  //! The minimum distance from the center to any edge of the bound.
  ElemType minimumBoundDistance;
  //! The dataset, consisting of end points
//  MatType* dataset;
  //! Collection of all the Segments in the tree
  std::vector<arma::mat>* dataset;

  //! The dataset consisting of Segments (their restrictions) located in this node.
  std::vector<arma::mat> & currentSegments

  //! Mapping between those two:
  std::vector<size_t> oldFromNew;
  //! M

 public:
  //! A single-tree traverser for binary space trees; see
  //! single_tree_traverser.hpp for implementation.
  template<typename RuleType>
  class SingleTreeTraverser;

  //! A dual-tree traverser for binary space trees; see dual_tree_traverser.hpp.
  template<typename RuleType>
  class DualTreeTraverser;

  template<typename RuleType>
  class BreadthFirstDualTreeTraverser;

  /**
   * Construct this as the root node of a binary space tree using the given
   * dataset.  This will copy the input matrix; if you don't want this, consider
   * using the constructor that takes an rvalue reference and use std::move().
   *
   * @param data Dataset of Segments to create tree from.  This will be copied!
   * @param maxLeafSize Size of each leaf in the tree.
   */
  BinarySpaceSegmentTree(const std::vector<arma::mat> & segments, const size_t maxLeafSize = 20);


  /**
   * Construct this node as a child of the given parent, starting at column
   * begin and using count points.  The ordering of that subset of points in the
   * parent's data matrix will be modified!  This is used for recursive
   * tree-building by the other constructors which don't specify point indices.
   *
   * A mapping of the old point indices to the new point indices is filled, but
   * it is expected that the vector is already allocated with size greater than
   * or equal to (begin + count), and if that is not true, invalid memory reads
   * (and writes) will occur.
   *
   * @param parent Parent of this node.  Its dataset will be modified!
   * @param oldFromNew Vector which will be filled with the old positions for
   *     each new point.
   * @param splitter Instantiated node splitter object.
   * @param maxLeafSize Size of each leaf in the tree.
   */
  BinarySpaceSegmentTree(BinarySpaceTree* parent,
                  std::vector<size_t>& oldFromNew,
                  SplitType<BoundType<MetricType>, MatType>& splitter,
                  const size_t maxLeafSize = 20);


  /**
   * Deletes this node, deallocating the memory for the children and calling
   * their destructors in turn.  This will invalidate any pointers or references
   * to any nodes which are children of this one.
   */
  ~BinarySpaceSegmentTree();

  //! Return the bound object for this node.
  const BoundType<MetricType>& Bound() const { return bound; }
  //! Return the bound object for this node.
  BoundType<MetricType>& Bound() { return bound; }

  //! Return the statistic object for this node.
  const StatisticType& Stat() const { return stat; }
  //! Return the statistic object for this node.
  StatisticType& Stat() { return stat; }

  //! Return whether or not this node is a leaf (true if it has no children).
  bool IsLeaf() const;

  //! Gets the left child of this node.
  BinarySpaceTree* Left() const { return left; }
  //! Modify the left child of this node.
  BinarySpaceTree*& Left() { return left; }

  //! Gets the right child of this node.
  BinarySpaceTree* Right() const { return right; }
  //! Modify the right child of this node.
  BinarySpaceTree*& Right() { return right; }

  //! Gets the parent of this node.
  BinarySpaceTree* Parent() const { return parent; }
  //! Modify the parent of this node.
  BinarySpaceTree*& Parent() { return parent; }

  //! Get the dataset which the tree is built on.
  const MatType& Dataset() const { return *dataset; }
  //! Modify the dataset which the tree is built on.  Be careful!
  MatType& Dataset() { return *dataset; }

  //! Get the metric that the tree uses.
  MetricType Metric() const { return MetricType(); }

  //! Return the number of children in this node.
  size_t NumChildren() const;

  /**
   * Return the index of the nearest child node to the given query point.  If
   * this is a leaf node, it will return NumChildren() (invalid index).
   */
  template<typename VecType>
  size_t GetNearestChild(
      const VecType& point,
      typename std::enable_if_t<IsVector<VecType>::value>* = 0);

  /**
   * Return the index of the furthest child node to the given query point.  If
   * this is a leaf node, it will return NumChildren() (invalid index).
   */
  template<typename VecType>
  size_t GetFurthestChild(
      const VecType& point,
      typename std::enable_if_t<IsVector<VecType>::value>* = 0);

  /**
   * Return the index of the nearest child node to the given query node.  If it
   * can't decide, it will return NumChildren() (invalid index).
   */
  size_t GetNearestChild(const BinarySpaceTree& queryNode);

  /**
   * Return the index of the furthest child node to the given query node.  If it
   * can't decide, it will return NumChildren() (invalid index).
   */
  size_t GetFurthestChild(const BinarySpaceTree& queryNode);

  /**
   * Return the furthest distance to a point held in this node.  If this is not
   * a leaf node, then the distance is 0 because the node holds no points.
   */
  ElemType FurthestPointDistance() const;

  /**
   * Return the furthest possible descendant distance.  This returns the maximum
   * distance from the centroid to the edge of the bound and not the empirical
   * quantity which is the actual furthest descendant distance.  So the actual
   * furthest descendant distance may be less than what this method returns (but
   * it will never be greater than this).
   */
  ElemType FurthestDescendantDistance() const;

  //! Return the minimum distance from the center of the node to any bound edge.
  ElemType MinimumBoundDistance() const;

  //! Return the distance from the center of this node to the center of the
  //! parent node.
  ElemType ParentDistance() const { return parentDistance; }
  //! Modify the distance from the center of this node to the center of the
  //! parent node.
  ElemType& ParentDistance() { return parentDistance; }

  /**
   * Return the specified child (0 will be left, 1 will be right).  If the index
   * is greater than 1, this will return the right child.
   *
   * @param child Index of child to return.
   */
  BinarySpaceTree& Child(const size_t child) const;

  BinarySpaceTree*& ChildPtr(const size_t child)
  { return (child == 0) ? left : right; }

  //! Return the number of points in this node (0 if not a leaf).
  size_t NumPoints() const;

  /**
   * Return the number of descendants of this node.  For a non-leaf in a binary
   * space tree, this is the number of points at the descendant leaves.  For a
   * leaf, this is the number of points in the leaf.
   */
  size_t NumDescendants() const;

  /**
   * Return the index (with reference to the dataset) of a particular descendant
   * of this node.  The index should be greater than zero but less than the
   * number of descendants.
   *
   * @param index Index of the descendant.
   */
  size_t Descendant(const size_t index) const;

  /**
   * Return the index (with reference to the dataset) of a particular point in
   * this node.  This will happily return invalid indices if the given index is
   * greater than the number of points in this node (obtained with NumPoints())
   * -- be careful.
   *
   * @param index Index of point for which a dataset index is wanted.
   */
  size_t Point(const size_t index) const;

  //! Return the minimum distance to another node.
  ElemType MinDistance(const BinarySpaceTree& other) const
  {
    return bound.MinDistance(other.Bound());
  }

  //! Return the maximum distance to another node.
  ElemType MaxDistance(const BinarySpaceTree& other) const
  {
    return bound.MaxDistance(other.Bound());
  }

  //! Return the minimum and maximum distance to another node.
  math::RangeType<ElemType> RangeDistance(const BinarySpaceTree& other) const
  {
    return bound.RangeDistance(other.Bound());
  }

  //! Return the minimum distance to another point.
  template<typename VecType>
  ElemType MinDistance(const VecType& point,
                       typename std::enable_if_t<IsVector<VecType>::value>* = 0)
      const
  {
    return bound.MinDistance(point);
  }

  //! Return the maximum distance to another point.
  template<typename VecType>
  ElemType MaxDistance(const VecType& point,
                       typename std::enable_if_t<IsVector<VecType>::value>* = 0)
      const
  {
    return bound.MaxDistance(point);
  }

  //! Return the minimum and maximum distance to another point.
  template<typename VecType>
  math::RangeType<ElemType>
  RangeDistance(const VecType& point,
                typename std::enable_if_t<IsVector<VecType>::value>* = 0) const
  {
    return bound.RangeDistance(point);
  }


  //! Store the center of the bounding region in the given vector.
  void Center(arma::vec& center) const { bound.Center(center); }

 private:
  /**
   * Splits the current node, assigning its left and right children recursively.
   *
   * @param maxLeafSize Maximum number of points held in a leaf.
   * @param splitter Instantiated SplitType object.
   */
  void SplitNode(const size_t maxLeafSize,
                 SplitType<BoundType<MetricType>, MatType>& splitter);

  /**
   * Update the bound of the current node. This method does not take into
   * account bound-specific properties.
   *
   * @param boundToUpdate The bound to update.
   * @param mt The end points of segments
   */
  template<typename BoundType2>
  void UpdateBound(BoundType2& boundToUpdate,arma::mat & mt);


};

} // namespace tree
} // namespace mlpack

// Include implementation.
#include "binary_space_segment_tree_impl.hpp"

// Include everything else, if necessary.
#include "../binary_space_tree.hpp"

#endif

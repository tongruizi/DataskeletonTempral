/**
 * @file binary_space_tree_impl.hpp
 *
 * Implementation of generalized space partitioning tree.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef MLPACK_CORE_TREE_BINARY_SPACE_SEGMENT_TREE_BINARY_SPACE_TREE_IMPL_HPP
#define MLPACK_CORE_TREE_BINARY_SPACE_SEGMENT_TREE_BINARY_SPACE_TREE_IMPL_HPP

// In case it wasn't included already for some reason.
#include "binary_space_segment_tree.hpp"

#include <mlpack/core/util/cli.hpp>
#include <mlpack/core/util/log.hpp>
#include <queue>

namespace mlpack
{
namespace tree
{

// Each of these overloads is kept as a separate function to keep the overhead
// from the two std::vectors out, if possible.

//! New constructor

template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType, SplitType>::
BinarySpaceSegmentTree(const std::vector<arma::mat> & SegmentData, const size_t maxLeafSize) :
    left(NULL),
    right(NULL),
    parent(NULL),
    bound(SegmentData[0].n_rows),
    parentDistance(0), // Parent distance for the root is 0: it has no parent.
    dataset(new std::vector<arma::mat>(SegmentData)), // Copies the dataset.
    currentSegments(SegmentData),
{
    //! Initilize oldFromNew
    (this->oldFromNew).resize(SegmentData.size());
    for (int i = 0; i < SegmentData.size(); i++)
    {
        (this->oldFromNew)[i] = i;
    }
    //! Do the actual splitting of this node.
    SplitType<BoundType<MetricType>, MatType> splitter;
    SplitNode(maxLeafSize, splitter);

    //! Create the statistic depending on if we are a leaf or not.
    stat = StatisticType(*this);
}


template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType, SplitType>::
BinarySpaceSegmentTree(
    BinarySpaceTree* parent,
    std::vector<arma::mat> & SegmentData,
    std::vector<size_t>& oldFromNew,
    SplitType<BoundType<MetricType>, MatType>& splitter,
    const size_t maxLeafSize) :
    left(NULL),
    right(NULL),
    parent(parent),
    bound(parent->Dataset().n_rows),
    dataset(&parent->Dataset()),
    currentSegments(SegmentData),
    oldFromNew(oldFromNew)
{
    // Hopefully the vector is initialized correctly!  We can't check that
    // entirely but we can do a minor sanity check.
    assert(oldFromNew.size() == dataset->n_cols);

    // Perform the actual splitting.
    SplitNode(oldFromNew, maxLeafSize, splitter);

    // Create the statistic depending on if we are a leaf or not.
    stat = StatisticType(*this);
}


/**
 * Deletes this node, deallocating the memory for the children and calling their
 * destructors in turn.  This will invalidate any pointers or references to any
 * nodes which are children of this one.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType, SplitType>::
~BinarySpaceSegmentTree()
{
    delete left;
    delete right;

    // If we're the root, delete the matrix.
    if (!parent)
        delete dataset;
}

template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline bool BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::IsLeaf() const
{
    return !left;
}

/**
 * Returns the number of children in this node.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::NumChildren() const
{
    if (left && right)
        return 2;
    if (left)
        return 1;

    return 0;
}

/**
 * Return the index of the nearest child node to the given query point.  If
 * this is a leaf node, it will return NumChildren() (invalid index).
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
template<typename VecType>
size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::GetNearestChild(
           const VecType& point,
           typename std::enable_if_t<IsVector<VecType>::value>*)
{
    if (IsLeaf() || !left || !right)
        return 0;

    if (left->MinDistance(point) <= right->MinDistance(point))
        return 0;
    return 1;
}

/**
 * Return the index of the furthest child node to the given query point.  If
 * this is a leaf node, it will return NumChildren() (invalid index).
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
template<typename VecType>
size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::GetFurthestChild(
           const VecType& point,
           typename std::enable_if_t<IsVector<VecType>::value>*)
{
    if (IsLeaf() || !left || !right)
        return 0;

    if (left->MaxDistance(point) > right->MaxDistance(point))
        return 0;
    return 1;
}

/**
 * Return the index of the nearest child node to the given query node.  If it
 * can't decide, it will return NumChildren() (invalid index).
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::GetNearestChild(const BinarySpaceTree& queryNode)
{
    if (IsLeaf() || !left || !right)
        return 0;

    ElemType leftDist = left->MinDistance(queryNode);
    ElemType rightDist = right->MinDistance(queryNode);
    if (leftDist < rightDist)
        return 0;
    if (rightDist < leftDist)
        return 1;
    return NumChildren();
}

/**
 * Return the index of the furthest child node to the given query node.  If it
 * can't decide, it will return NumChildren() (invalid index).
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::GetFurthestChild(const BinarySpaceTree& queryNode)
{
    if (IsLeaf() || !left || !right)
        return 0;

    ElemType leftDist = left->MaxDistance(queryNode);
    ElemType rightDist = right->MaxDistance(queryNode);
    if (leftDist > rightDist)
        return 0;
    if (rightDist > leftDist)
        return 1;
    return NumChildren();
}

/**
 * Return a bound on the furthest point in the node from the center.  This
 * returns 0 unless the node is a leaf.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline
typename BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::ElemType
         BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::FurthestPointDistance() const
{
    if (!IsLeaf())
        return 0.0;

    // Otherwise return the distance from the center to a corner of the bound.
    return 0.5 * bound.Diameter();
}

/**
 * Return the furthest possible descendant distance.  This returns the maximum
 * distance from the center to the edge of the bound and not the empirical
 * quantity which is the actual furthest descendant distance.  So the actual
 * furthest descendant distance may be less than what this method returns (but
 * it will never be greater than this).
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline
typename BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::ElemType
         BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::FurthestDescendantDistance() const
{
    return furthestDescendantDistance;
}

//! Return the minimum distance from the center to any bound edge.
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline
typename BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::ElemType
         BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
         SplitType>::MinimumBoundDistance() const
{
    return bound.MinWidth() / 2.0;
}

/**
 * Return the specified child.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>&
       BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::Child(const size_t child) const
{
    if (child == 0)
        return *left;
    else
        return *right;
}

/**
 * Return the number of points contained in this node.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::NumPoints() const
{
    if (left)
        return 0;

    return this->currentSegments.size();
}

/**
 * Return the number of descendants contained in the node.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::NumDescendants() const
{
    return this->currentSegments.size();
}

/**
 * Return the index of a particular descendant contained in this node.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::Descendant(const size_t index) const
{
    return this->oldFromNew[index];
}

/**
 * Return the index of a particular point contained in this node.
 */
template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
inline size_t BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType,
       SplitType>::Point(const size_t index) const
{
    return this->oldFromNew[index];
}


template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
void BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType, SplitType>::
SplitNode(const size_t maxLeafSize,
          SplitType<BoundType<MetricType>, MatType>& splitter)
{
    // We need to expand the bounds of this node properly.
    arma::mat boundaryPoints;
    PointExtractor::ExtractPoints(boundaryPoints, this->currentSegments);

    UpdateBound(bound);

    // Calculate the furthest descendant distance.
    furthestDescendantDistance = 0.5 * bound.Diameter();

    // First, check if we need to split at all.
    if (this->currentSegments.size() <= maxLeafSize)
        return; // We can't split this.

    // splitCol denotes the two partitions of the dataset after the split. The
    // points on its left go to the left child and the others go to the right
    // child.
//    size_t splitCol;

    // Find the partition of the node. This method does not perform the split.
    typename Split::SplitInfo splitInfo;

    const bool split = splitter.SplitNode(bound, boundaryPoints, splitInfo);

    // The node may not be always split. For instance, if all the points are the
    // same, we can't split them.
    if (!split)
        return;

    // Perform the actual splitting.  This will order the dataset such that
    // points that belong to the left subtree are on the left of splitCol, and
    // points from the right subtree are on the right side of splitCol.

    std::vector<arma::mat> segmentLeft;
    std::vector<arma::mat> segmentRight;
    std::vector<size_t> oldFromNewLeft;
    std::vector<size_t> oldFromNewRight;

    splitter.PerformSplit(this->currentSegments, segmentLeft, segmentRight,
    splitInfo, this->oldFromNew, oldFromNewLeft, oldFromNewRight);

   // assert(splitCol > begin);
  //  assert(splitCol < begin + count);

    // Now that we know the split column, we will recursively split the children
    // by calling their constructors (which perform this splitting process).
    left = new BinarySpaceTree(this, segmentLeft, oldFromNewLeft, splitter, maxLeafSize);
    right = new BinarySpaceTree(this,segmentRight, oldFromNewRight, splitter, maxLeafSize);

    // Calculate parent distances for those two nodes.
    arma::vec center, leftCenter, rightCenter;
    Center(center);
    left->Center(leftCenter);
    right->Center(rightCenter);

    const ElemType leftParentDistance = MetricType::Evaluate(center, leftCenter);
    const ElemType rightParentDistance = MetricType::Evaluate(center, rightCenter);

    left->ParentDistance() = leftParentDistance;
    right->ParentDistance() = rightParentDistance;
}

template<typename MetricType,
         typename StatisticType,
         typename MatType,
         template<typename BoundMetricType, typename...> class BoundType,
         template<typename SplitBoundType, typename SplitMatType>
         class SplitType>
template<typename BoundType2>
void BinarySpaceSegmentTree<MetricType, StatisticType, MatType, BoundType, SplitType>::
UpdateBound(BoundType2& boundToUpdate,arma::mat & mt)
{
        boundToUpdate |= mt;
}



} // namespace tree
} // namespace mlpack

#endif

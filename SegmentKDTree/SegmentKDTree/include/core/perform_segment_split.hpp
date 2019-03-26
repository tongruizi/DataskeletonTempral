/**
 * @file perform_segment_split.hpp
 * @author Yury Elkin
 *
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef MLPACK_CORE_TREE_PERFORM_SEGMENT_SPLIT_HPP
#define MLPACK_CORE_TREE_PERFORM_SEGMENT_SPLIT_HPP

namespace mlpack
{
namespace tree /** Trees and tree-building procedures. */
{
namespace split
{

/**
 * This function implements the default split behavior i.e. it rearranges
 * points according to the split information. The SplitType::AssignToLeftNode()
 * function is used in order to determine the child that contains any particular
 * point.
 *
 * @param data The dataset used by the binary space tree.
 * @param begin Index of the starting point in the dataset that belongs to
 *    this node.
 * @param count Number of points in this node.
 * @param splitInfo The information about the split.
 */
template<typename MatType, typename SplitType>
void PerformSegmentSplit(std::vector<arma::mat> & input, std::vector<arma::mat> & leftout,
                  std::vector<arma::mat> & rightout,
                  const typename SplitType::SplitInfo& splitInfo)
{
    // This method modifies the input dataset.  We loop both from the left and
    // right sides of the points contained in this node.

    int dim = input[i].n_rows;
    for (int i = 0; i < intput.size; i++)
    {
        bool f1 = SplitType::AssignToLeftNode(input[i].col(0), splitInfo);
        bool f2 = SplitType::AssignToLeftNode(input[i].col(1), splitInfo);
        if (f1 && f2)
        {
            leftout.push_back(input[i]);
        }
        else if ((!f1) && (!f2))
        {
            rightout.push_back(input[i]);
        }
        else if (f1 && (!f2))
        {
            arma::mat leftSegment(dim,2);
            arma::mat rightSegment(dim,2);
            arma::vec newPoint;
            SplitType::FindSegmentSplit(input[i].col(0),input[i].col(1), newPoint,
                                        splitInfo);
            leftSegment.col(0) = input[i].col(0);
            leftSegment.col(1) = newPoint;
            rightSegment.col(0) = newPoint;
            rightSegment.col(1) = input[i].col(1);
            leftout.push_back(leftSegment);
            rightout.push_back(rightSegment);

        }
        else if ((!f1) && f2)
        {
            arma::mat leftSegment(dim,2);
            arma::mat rightSegment(dim,2);
            arma::vec newPoint;
            SplitType::FindSegmentSplit(input[i].col(0),input[i].col(1), newPoint,
                                        splitInfo);
            leftSegment.col(0) = input[i].col(1);
            leftSegment.col(1) = newPoint;
            rightSegment.col(0) = newPoint;
            rightSegment.col(1) = input[i].col(0);
            leftout.push_back(leftSegment);
            rightout.push_back(rightSegment);

        }


    }

}

/**
 * This function implements the default split behavior i.e. it rearranges
 * points according to the split information. The SplitType::AssignToLeftNode()
 * function is used in order to determine the child that contains any particular
 * point. The function takes care of indices and returns the list of changed
 * indices.
 *
 * @param data The dataset used by the binary space tree.
 * @param begin Index of the starting point in the dataset that belongs to
 *    this node.
 * @param count Number of points in this node.
 * @param splitInfo The information about the split.
 * @param oldFromNew Vector which will be filled with the old positions for
 *    each new point.
 */
template<typename MatType, typename SplitType>
size_t PerformSegmentSplit(std::vector<arma::mat> & input, std::vector<arma::mat> & leftout,
                    std::vector<arma::mat> & rightout,
                    const typename SplitType::SplitInfo& splitInfo,
                    std::vector<size_t>& oldFromNew
                    std::vector<size_t>& oldFromNewLeft,
                    std::vector<size_t>& oldFromNewRight)
{
    // This method modifies the input dataset.  We loop both from the left and
    // right sides of the points contained in this node.

    int dim = input[i].n_rows;
    for (int i = 0; i < intput.size; i++)
    {
        bool f1 = SplitType::AssignToLeftNode(input[i].col(0), splitInfo);
        bool f2 = SplitType::AssignToLeftNode(input[i].col(1), splitInfo);
        if (f1 && f2)
        {
            leftout.push_back(input[i]);
            oldFromNewLeft.push_back(oldFromNew[i]);
        }
        else if ((!f1) && (!f2))
        {
            rightout.push_back(input[i]);
        }
        else if (f1 && (!f2))
        {
            arma::mat leftSegment(dim,2);
            arma::mat rightSegment(dim,2);
            arma::vec newPoint;
            SplitType::FindSegmentSplit(input[i].col(0),input[i].col(1), newPoint,
                                        splitInfo);
            leftSegment.col(0) = input[i].col(0);
            leftSegment.col(1) = newPoint;
            rightSegment.col(0) = newPoint;
            rightSegment.col(1) = input[i].col(1);
            leftout.push_back(leftSegment);
            oldFromNewLeft.push_back(oldFromNew[i]);
            rightout.push_back(rightSegment);
            oldFromNewRight.push_back(oldFromNew[i]);

        }
        else if ((!f1) && f2)
        {
            arma::mat leftSegment(dim,2);
            arma::mat rightSegment(dim,2);
            arma::vec newPoint;
            SplitType::FindSegmentSplit(input[i].col(0),input[i].col(1), newPoint,
                                        splitInfo);
            leftSegment.col(0) = input[i].col(1);
            leftSegment.col(1) = newPoint;
            rightSegment.col(0) = newPoint;
            rightSegment.col(1) = input[i].col(0);
            leftout.push_back(leftSegment);
            oldFromNewLeft.push_back(oldFromNew[i]);
            rightout.push_back(rightSegment);
            oldFromNewRight.push_back(oldFromNew[i]);
        }

    }
}


} // namespace split
} // namespace tree
} // namespace mlpack


#endif // MLPACK_CORE_TREE_BINARY_SPACE_TREE_PERFORM_SPLIT_HPP

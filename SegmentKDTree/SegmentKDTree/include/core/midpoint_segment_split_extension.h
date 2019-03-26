//#ifndef MIDPOINT_SEGMENT_SPLIT_EXTENSION_H
//#define MIDPOINT_SEGMENT_SPLIT_EXTENSION_H
//
//#include <mlpack/core/tree/binary_space_tree/midpoint_split.hpp>
//template<typename BoundType, typename MatType = arma::mat>
//class midpoint_segment_split_extension : public mlpack::tree::MidpointSplit<BoundType,MatType>
//{
//public:
//    midpoint_segment_split_extension() {}
//    static void FindSegmentSplit(const VecType& pointA, const VecType& pointB, const VecType & result
//                     const SplitInfo& splitInfo)
//    {
//        result.set_size(pointA.n_rows,pointA.n_cols);
//        int dimension = splitInfo.splitDimension;
//        double t = (splitInfo.splitVal - pointA(dimension)) / (pointB(dimension) - pointA(dimension));
//        result = pointA + t*(pointA - pointB);
//    }
//
//
//protected:
//
//private:
//};
//
//#endif // MIDPOINT_SEGMENT_SPLIT_EXTENSINO_H

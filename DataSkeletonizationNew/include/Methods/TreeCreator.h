#ifndef TREECREATOR_H
#define TREECREATOR_H


class TreeCreator
{
public:
    TreeCreator() {}
    template<typename TreeType, typename MatType>
    static TreeType* BuildTree(
        MatType&& dataset,
        std::vector<size_t>& oldFromNew,
        const typename std::enable_if<
        mlpack::tree::TreeTraits<TreeType>::RearrangesDataset>::type* = 0)
    {
        return new TreeType(std::forward<MatType>(dataset), oldFromNew);
    }

//! Call the tree constructor that does not do mapping.
    template<typename TreeType, typename MatType>
    static TreeType* BuildTree(
        MatType&& dataset,
        const std::vector<size_t>& /* oldFromNew */,
        const typename std::enable_if<
        !mlpack::tree::TreeTraits<TreeType>::RearrangesDataset>::type* = 0)
    {
        return new TreeType(std::forward<MatType>(dataset));
    }

protected:

private:
};

#endif // TREECREATOR_H

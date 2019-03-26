#ifndef TREEEXPEREMENTS_H
#define TREEEXPEREMENTS_H

#include <iostream>
#include <mlpack/core.hpp>
#include <mlpack/prereqs.hpp>
#include <mlpack/core/tree/cover_tree.hpp>

class TreeExperements
{
    public:
        TreeExperements() {}
        void Test_Cover_Tree(arma::mat & data)
        {
        mlpack::tree::CoverTree<> treep(data);
        std::cout << "Descendant number: " << treep.NumDescendants() << std::endl;
        }
        void Test_KD_Tree(arma::mat & data)
        {
  //      mlpack::tree::BinarySpaceTree<metric::EuclideanDistance>(data);

        }

    protected:

    private:
};

#endif // TREEEXPEREMENTS_H

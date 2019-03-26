#ifndef TOSPTTRANSFORMER_H
#define TOSPTTRANSFORMER_H

#include <mlpack/core.hpp>

struct VectorComparator
{
    bool operator() (arma::Col<int>  a, arma::Col<int>  b) const
    {
        for (int i = 0; i < a.size(); i++)
        {
            if (a(i) > b(i))
            {
                return true;
            }
            else if (b(i) < a(i))
            {
                return false;
            }

        }
        return true;
    }
};

class ToSptTransformer
{
public:
    ToSptTransformer() {}

    static int minus2odd_or_even(int i)
    {
        if (i % 2 == 0)
        {
            return 1;

        }

        else
        {
            return -1;
        }


    }

    static void buildTetraGrid(int nx, int ny, int nz, arma::Mat<int> & tetra,arma::Cube<int> & dicIJK2Index)
    {
        int STATE1 = 1;
        int STATE2 = -1;
        int nTetra = (nx-1) * (ny-1) * (nz -1) * 5;
        tetra.set_size(4,nTetra);
        int index_tetra = 0;
        int cur_state = STATE1;

        for (int i = 0; i < nx-1; i++)
        {
            for (int j = 0 ; j < ny-1; j++)
            {
                cur_state = STATE1 * minus2odd_or_even(i + j);
                for (int k  = 0; k < nz-1; k++)
                {
                    int pt1_ind = dicIJK2Index(i, j, k);
                    int pt2_ind = dicIJK2Index(i, j + 1, k);
                    int pt3_ind = dicIJK2Index(i + 1, j, k);
                    int pt4_ind = dicIJK2Index(i + 1, j + 1, k);
                    int pt5_ind = dicIJK2Index(i, j, k + 1);
                    int pt6_ind = dicIJK2Index(i, j + 1, k + 1);
                    int pt7_ind = dicIJK2Index(i + 1, j, k + 1);
                    int pt8_ind = dicIJK2Index(i + 1, j + 1, k + 1);
                    if (cur_state == STATE1)
                    {
                        arma::Col<int> w{pt1_ind, pt2_ind, pt3_ind, pt5_ind};
                        tetra.col(index_tetra) = w;
                        arma::Col<int> w2 {pt3_ind, pt5_ind, pt7_ind, pt8_ind};
                        tetra.col(index_tetra+1) = w2;
                        arma::Col<int> w3 {pt2_ind, pt3_ind, pt4_ind, pt8_ind};
                        tetra.col(index_tetra+2) = w3;
                        arma::Col<int> w4 {pt2_ind, pt5_ind, pt6_ind, pt8_ind};
                        tetra.col(index_tetra+3) = w4;
                        arma::Col<int> w5 {pt2_ind, pt3_ind, pt5_ind, pt8_ind};
                        tetra.col(index_tetra+4) = w5;


                    }

                    if (cur_state == STATE2)
                    {
                        arma::Col<int> w {pt1_ind, pt3_ind, pt4_ind, pt7_ind};
                        tetra.col(index_tetra) = w;
                        arma::Col<int> w2 {pt1_ind, pt2_ind, pt4_ind, pt6_ind};
                        tetra.col(index_tetra+1) = w2;
                        arma::Col<int> w3 {pt4_ind, pt6_ind, pt7_ind, pt8_ind};
                        tetra.col(index_tetra+2) = w3;
                        arma::Col<int> w4 {pt1_ind, pt5_ind, pt6_ind, pt7_ind};
                        tetra.col(index_tetra+3) = w4;
                        arma::Col<int> w5 {pt1_ind, pt4_ind, pt6_ind, pt7_ind};
                        tetra.col(index_tetra+4) = w5;

                    }
                    index_tetra = index_tetra + 5;
                    cur_state = -cur_state;
                }
            }

        }
    }
    static void buildTriFromTetra(arma::Mat<int> & tetra, arma::Mat<int> & tri_array)
    {
        std::cout << "Tri transformation 1 " << std::endl;
        std::set<std::array<int,3>> tri;
        //std::set<arma::Col<int>,VectorComparator> tri;

        int nTe = tetra.n_cols;

        std::cout << "Tetra col size:" << nTe << std::endl;

        int tri_index = 0;
        for (int i = 0; i < nTe; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                //arma::Col<int> newTri(3);
                std::array<int,3> newTri;
                int li = 0;
                for (int k = 0; k < 4; k++)
                {
                    if (k != j)
                    {
                        newTri[li] = tetra(k,i);
                        li++;
                    }

                }
                std::sort(newTri.begin(), newTri.end());
                if (tri.find(newTri) == tri.end())
                {
                    tri.insert(newTri);
                }

            }


        }

        std::cout << "Tri transformation 2 " << std::endl;

        tri_array.set_size(3,tri.size());
        int indx = 0;
        std::cout << "Tri size right now : " << tri.size() << std::endl;
        for (auto it = tri.begin(); it != tri.end(); it++)
        {
            arma::Col<int> w(3);
            for (int k = 0; k < 3; k++)
            {
                w(k) = (*it)[k];
            }
            tri_array.col(indx) = w;
            indx++;
        }
    }
    static void buildEdgeFromTri(arma::Mat<int> & tri, arma::Mat<int> & edge_array)
    {
        std::cout << "Edgesize: " << tri.n_cols << std::endl;
        //    std::set<arma::Col<int>,VectorComparator> edge;
        std::set<std::array<int,2>> edge;
        int nTri = tri.n_cols;

        for (int i = 0; i < nTri; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                //arma::Col<int> newEdge(2);
                std::array<int,2> newEdge;
                int li = 0;
                for (int k = 0; k < 3; k++)
                {
                    if (k != j)
                    {
                        newEdge[li] = tri(k,i);
                        li++;
                    }

                }
                std::sort(newEdge.begin(), newEdge.end());
                if (edge.find(newEdge) == edge.end())
                {
                    edge.insert(newEdge);
                }

            }
        }
        edge_array.set_size(2,edge.size());
        int indx = 0;
        std::cout << "Edge totalsize: " << edge.size() << std::endl;
        for (auto it = edge.begin(); it != edge.end(); it++)
        {
            arma::Col<int> w(2);
            for (int k = 0; k < 2; k++)
            {
                w(k) = (*it)[k];
            }
            edge_array.col(indx) = w;
            indx++;
        }
    }
    // /home/yury/LocalTests/DebugDebug
    static void PrintMatrixToFile(std::string & w, arma::Mat<int> matrix)
    {
        std::ofstream mystream;
        mystream.open(w);
        matrix = matrix.t();
        mystream << matrix << std::endl;
       // for (int i = 0; i < matrix.n_cols; i++)
      //  {
      //      for (int j = 0; j < matrix.n_rows; j++)
     //       {
        //        mystream << matrix(j,i) << std::endl;
      //      }
     //   }
        mystream.close();

    }




protected:

private:
};

#endif // TOSPTTRANSFORMER_H

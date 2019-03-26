#ifndef DIMOSCLAUNCHER_H
#define DIMOSCLAUNCHER_H

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
using namespace std;

#include "Simplex.h"
#include "persistence.h"
#include "DiscreteVField.h"
#include "Simplicial2Complex.h"
#include "Graph.h"


class DiMoSCLauncher
{
public:
    DiMoSCLauncher() {}
    static void Launch(MyGraphType & G, double ve_delta, arma::mat & data, std::vector<double> & f, arma::Mat<int> & edges, arma::Mat<int> & triangles)
    {
        clock_t begin = clock();
        int DIM = data.n_rows;
        //------------------- pre_time ------------------------
        string pre_save;
        bool use_pre_save = false;

       // double et_delta = 0;
      //  double ve_delta = 0;

        //    void buildComplexFromInfo(arma::mat & points, std::vector<double> & f, arma::Mat<int> & tri, arma::Mat<int> & edges);

      //  ve_delta = atof(argv[3]);
      //  DIM = atoi(argv[4]);


        Simplicial2Complex K;
        clock_t begin_preprocess = clock();
        if (!use_pre_save)
        {

            clock_t begin_ = clock();

           //cout << "Reading in simplicial complex...\n";
            //K.buildComplexFromFile2_BIN(argv[1]);
            K.buildComplexFromInfo(data,f,triangles,edges);
           // cout << "Done\n";
           // cout.flush();
            // cin.get(); // has test

            clock_t end_ = clock();
            double elapsed_secs_ = double(end_ - begin_) / CLOCKS_PER_SEC;
           // cout<<"read in and build complex: "<<elapsed_secs_<<endl;

            begin_ = clock();
            // Build psudo morse function
           // cout << "Building pseudo-Morse function...\n";
            K.buildPsuedoMorseFunction();
          //  cout << "Done\n";
          //  cout.flush();
            // cin.get();

            end_ = clock();

            elapsed_secs_ = double(end_ - begin_) / CLOCKS_PER_SEC;
         //   cout<<"building pseudo-Morse function: "<<elapsed_secs_<<endl;

            begin_ = clock();

          //  cout << "Building filtration...\n";
            K.buildFiltrationWithLowerStar();
         //   cout << "Done\n";
         //   cout.flush();
            // cin.get(); // has test

            end_ = clock();
            elapsed_secs_ = double(end_ - begin_) / CLOCKS_PER_SEC;
         //   cout<<"filtration time: "<<elapsed_secs_<<endl;

            begin_ = clock();

          //  cout << "Computing persistence pairs...\n";
            K.PhatPersistence();
          //  cout << "Done!\n";
          //  cout.flush();
          //  cout << "Outputing persistence pairs...\n";

            end_ = clock();
            elapsed_secs_ = double(end_ - begin_) / CLOCKS_PER_SEC;
           // cout<<"persistence time: "<<elapsed_secs_<<endl;

           // cout << "Writing pre_saved_data...\n";
         //   K.write_presave(argv[2]);
           // cout << "Done!\n";
           // cout.flush();
        }
        else
        {
           // cout << "Reading in pre_saved_data...\n";
           // K.Load_Presaved(argv[1], pre_save);
           // cout << "Done\n";
           // cout.flush();
            // cin.get();

            // Build psudo morse function
           // cout << "Building pseudo-Morse function...\n";
            K.buildPsuedoMorseFunction();
          //  cout << "Done\n";
          //  cout.flush();
        }
        //------------------- pre_time ------------------------
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

        clock_t end_preprocess = clock();
        double preprocess_time = end_preprocess-begin_preprocess;
        // cout<<"preprocess time: "<<preprocess_time/ CLOCKS_PER_SEC<<endl;

//sp2
        clock_t begin3 = clock();
// build tree
        // cout<<"sp2:"<<endl;
        begin = clock();

       // cout<<"build spanning tree.";
        K.build_spt3(ve_delta);
       // cout<<"Done."<<endl;

        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       // cout<<"build sp tree time: "<<elapsed_secs<<endl;
//retrive
        begin = clock();

        //K.retrieve_1manifold2();
        K.retrieve_1manifold_simp();

        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       // cout<<"retrive time: "<<elapsed_secs<<endl;
//output
        begin = clock();

        K.convert_output();
        K.write_to_graph(G);

        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      //  cout<<"write output time: "<<elapsed_secs<<endl;
//sp2
        clock_t end3 = clock();
        double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
       // cout<<"sp2 time: "<<elapsed_secs3<<endl;



    }
    virtual ~DiMoSCLauncher() {}

protected:

private:
};

#endif // DIMOSCLAUNCHER_H

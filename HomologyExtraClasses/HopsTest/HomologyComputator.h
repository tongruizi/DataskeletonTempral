#ifndef HOMOLOGYCOMPUTATOR_H
#define HOMOLOGYCOMPUTATOR_H

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <iostream>
#include <ctime>
#include <utility>
#include <vector>



class HomologyComputator
{
    using Simplex_tree = Gudhi::Simplex_tree<>;
    using Filtration_value = Simplex_tree::Filtration_value;
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;
    using typeVectorVertex = std::vector< Simplex_tree::Vertex_handle >;

public:


    HomologyComputator() {}

    static void ComputeOriginal(Simplex_tree & st, int coeff_field_characteristic)
    {
        Filtration_value min_persistence = 0.0;
        Persistent_cohomology pcoh(st);
        // initializes the coefficient field for homology
        pcoh.init_coefficients(coeff_field_characteristic);
        pcoh.compute_persistent_cohomology(min_persistence);
        // Output the diagram in filediag
        pcoh.output_diagram();

    }

    static void Compute(int coeff_field_characteristic)
    {
        // program args management

        Filtration_value min_persistence = 0.0;

        // TEST OF INSERTION
        std::cout << "********************************************************************" << std::endl;
        std::cout << "TEST OF INSERTION" << std::endl;
        Simplex_tree st;
        // ++ FIRST
        std::cout << "   - INSERT (0,1,2)" << std::endl;
        typeVectorVertex SimplexVector = {0, 1, 2};
        st.insert_simplex_and_subfaces(SimplexVector, 0.3);
        // ++ SECOND
        std::cout << "   - INSERT 3" << std::endl;
        SimplexVector = {3};
        st.insert_simplex_and_subfaces(SimplexVector, 0.1);
        // ++ THIRD
        std::cout << "   - INSERT (0,3)" << std::endl;
        SimplexVector = {0, 3};
        st.insert_simplex_and_subfaces(SimplexVector, 0.2);
        // ++ FOURTH
        std::cout << "   - INSERT (0,1) (already inserted)" << std::endl;
        SimplexVector = {0, 1};
        st.insert_simplex_and_subfaces(SimplexVector, 0.2);
        // ++ FIFTH
        std::cout << "   - INSERT (3,4,5)" << std::endl;
        SimplexVector = {3, 4, 5};
        st.insert_simplex_and_subfaces(SimplexVector, 0.3);
        // ++ SIXTH
        std::cout << "   - INSERT (0,1,6,7)" << std::endl;
        SimplexVector = {0, 1, 6, 7};
        st.insert_simplex_and_subfaces(SimplexVector, 0.4);
        // ++ SEVENTH
        std::cout << "   - INSERT (4,5,8,9)" << std::endl;
        SimplexVector = {4, 5, 8, 9};
        st.insert_simplex_and_subfaces(SimplexVector, 0.4);
        // ++ EIGHTH
        std::cout << "   - INSERT (9,10,11)" << std::endl;
        SimplexVector = {9, 10, 11};
        st.insert_simplex_and_subfaces(SimplexVector, 0.3);
        // ++ NINETH
        std::cout << "   - INSERT (2,10,12)" << std::endl;
        SimplexVector = {2, 10, 12};
        st.insert_simplex_and_subfaces(SimplexVector, 0.3);
        // ++ TENTH
        std::cout << "   - INSERT (11,6)" << std::endl;
        SimplexVector = {6, 11};
        st.insert_simplex_and_subfaces(SimplexVector, 0.2);
        // ++ ELEVENTH
        std::cout << "   - INSERT (13,14,15)" << std::endl;
        SimplexVector = {13, 14, 15};
        st.insert_simplex_and_subfaces(SimplexVector, 0.25);

        std::cout << "The complex contains " << st.num_simplices() << " simplices - " << st.num_vertices() << " vertices "
                  << std::endl;
        std::cout << "   - dimension " << st.dimension() << std::endl;
        std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:"
                  << std::endl;
        std::cout << "**************************************************************" << std::endl;
        std::cout << "strict graph G { " << std::endl;
        for (auto f_simplex : st.filtration_simplex_range())
        {
            std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
            for (auto vertex : st.simplex_vertex_range(f_simplex))
            {
                std::cout << static_cast<int>(vertex) << " -- ";
            }
            std::cout << ";" << std::endl;
        }
        std::cout << "}" << std::endl;
        std::cout << "**************************************************************" << std::endl;
        // Compute the persistence diagram of the complex
        Persistent_cohomology pcoh(st);
        // initializes the coefficient field for homology
        pcoh.init_coefficients(coeff_field_characteristic);
        pcoh.compute_persistent_cohomology(min_persistence);
        // Output the diagram in filediag
        pcoh.output_diagram();

    }
protected:

private:
};

#endif // HOMOLOGYCOMPUTATOR_H

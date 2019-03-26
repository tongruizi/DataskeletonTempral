#include "HomologyComputator.h"

void testSimplexMassiveAddition(Simplex_tree & st, int k)
{
    for (int i = 0; i < k; i++)
    {
        auto apm = {i};
        st.insert_simplex_and_subfaces(apm,i);
    }
    for (int i = 1; i < k; i++)
    {
        auto apm = {i,0};
        st.insert_simplex_and_subfaces(apm,2*i);
    }

}

HomologyComputator::HomologyComputator()
{
    //ctor
}

void HomologyComputator::TestRun()
{
    int coeff_field_characteristic = 2;
    int min_persistence = 0;
    Simplex_tree st;
    testSimplexMassiveAddition(st,10000);
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

void HomologyComputator::Compute(Simplex_tree & st)
{
    int coeff_field_characteristic = 2;
    int min_persistence = 0;
    Persistent_cohomology pcoh(st);
    pcoh.init_coefficients(coeff_field_characteristic);
    pcoh.compute_persistent_cohomology(min_persistence);
    pcoh.output_diagram();

}


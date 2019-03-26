#include "binomial_coeff_table.h"


binomial_coeff_table::binomial_coeff_table(index_t n, index_t k) {
		binomial_coeff_table::n_max = n;
		binomial_coeff_table::k_max = k;

		B.resize(n + 1);
		for (index_t i = 0; i <= n; i++) {
			B[i].resize(k + 1);
			for (index_t j = 0; j <= std::min(i, k); j++) {
				if (j == 0 || j == i)
					B[i][j] = 1;
				else
					B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			}
		}
	}

	index_t binomial_coeff_table::operator()(index_t n, index_t k) const {
    if (k > n)
    {
        return 0;
    }
    return B[n][k];
}

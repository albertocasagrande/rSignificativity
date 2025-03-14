#include <vector>
#include <cmath>

#include <gmpxx.h>

#include <Rcpp.h>

template <typename ELEM_TYPE>
using matrix = std::vector<std::vector<ELEM_TYPE>>;

/**
 * @brief Compute the binomial coefficient
 *
 * This function compute the binomial coefficient \$\binom{a}{b}\$.
 * According to the uniform cost criterion, its asymptotic time
 * complexity is \$O(\min\{b, a-b\})\$.
 *
 * @tparam PARAM_TYPE is the type of the parameters
 * @tparam OUTPUT_TYPE is the type of the output
 * @param a is an integer value
 * @param b is an integer value
 * @return OUTPUT_TYPE
 */
template <typename PARAM_TYPE, typename OUTPUT_TYPE = PARAM_TYPE>
OUTPUT_TYPE binom(PARAM_TYPE a, PARAM_TYPE b)
{
    if (a < b)
    {
        throw std::domain_error("The first parameter must be greater than the second one.");
    }

    PARAM_TYPE c{a - b};
    if (c > b)
    {
        c = b;
    }

    OUTPUT_TYPE result{1};
    for (PARAM_TYPE i = 0; i < c; ++i)
    {
        result = (result * (a - i)) / (i + 1);
    }

    return result;
}

/**
 * @brief A lexicographic order enumerator for weak compositions of m in k parts
 *
 * This function implements a lexicographic order enumerator for weak
 * compositions of m in k parts. A weak composition of m in k parts is a vector
 * of natural numbers whose sum is m. The lexicographic order \$<_l\$ of weak
 * compositions of m in k parts is an order such that
 * \$\langle a_1, \ldots a_k \rangle <_l \langle b_1, \ldots b_k \rangle\$ when
 * there exist \$i \in [1, k]\$ such that for all \$j \in [1, i-1]\$ it holds
 * that $a_j = b_j$ and $a_i<b_i$.
 *
 * This method returns the `i`-th weak composition in the lexicographic order
 * of weak compositions of m in k parts.
 *
 * @tparam VECTOR_TYPE is the weak component type
 * @tparam INDEX_TYPE is the index type
 * @param ell is the output vector
 * @param m is the sum of the output vector components
 * @param k is the size of the output vector
 * @param i is the index
 */
template <typename VECTOR_TYPE, typename INDEX_TYPE>
void fill_iota(VECTOR_TYPE &ell, unsigned int m, unsigned int k, INDEX_TYPE i)
{
    if (ell.size() < k)
    {
        throw std::domain_error("The ell must have k components at least");
    }
    auto ell_it = ell.begin();

    INDEX_TYPE C = binom<unsigned int, INDEX_TYPE>(m + k - 2, m);
    while (m > 0 && k > 1)
    {
        unsigned int ml = m;
        bool loop_cond = (i >= C);
        while (loop_cond)
        {
            i = i-C;

            loop_cond = false;
            if (m > 0)
            {
                INDEX_TYPE B = (m * C) / (m + k - 2);

                if (i >= B)
                {
                    C = B;
                    --m;
                    loop_cond = true;
                }
            }
        }

        *ell_it = ml - m;
        ++ell_it;

        --k;
        C = (k - 1) * C / (m + k - 1);
    }
    *ell_it = m;
}

/**
 * @brief A lexicographic order enumerator for weak compositions of m in k parts
 *
 * This function implements a lexicographic order enumerator for weak
 * compositions of m in k parts. A weak composition of m in k parts is a vector
 * of natural numbers whose sum is m. The lexicographic order \$<_l\$ of weak
 * compositions of m in k parts is an order such that
 * \$\langle a_1, \ldots a_k \rangle <_l \langle b_1, \ldots b_k \rangle\$ when
 * there exist \$i \in [1, k]\$ such that for all \$j \in [1, i-1]\$ it holds
 * that $a_j = b_j$ and $a_i<b_i$.
 *
 * This method returns the `i`-th weak composition in the lexicographic order
 * of weak compositions of m in k parts.
 *
 * @tparam VECTOR_TYPE is the weak component type
 * @tparam INDEX_TYPE is the index type
 * @param m is the sum of the output vector components
 * @param k is the size of the output vector
 * @param i is the index
 * @return std::vector<unsigned int>
 */
template <typename VECTOR_TYPE, typename INDEX_TYPE>
VECTOR_TYPE
iota(const unsigned int m, const unsigned int k, const INDEX_TYPE i)
{
    VECTOR_TYPE ell(k);

    fill_iota(ell, m, k, i);

    return ell;
}

/**
 * @brief Fill a square matrix with the components of a vector
 *
 * @tparam ELEM_TYPE is the type of the vector components
 * @param V is the input vector
 * @return matrix<ELEM_TYPE>
 */
template <typename ELEM_TYPE>
void fill_gamma(matrix<ELEM_TYPE> &M, const std::vector<ELEM_TYPE> &V)
{
    if (M.size() * M.size() != V.size())
    {
        throw std::domain_error("The vector cannot be reshaped in the square matrix");
    }

    auto V_it = V.begin();
    for (auto row_it = M.begin(); row_it != M.end(); ++row_it)
    {
        std::vector<ELEM_TYPE> &row = *row_it;
        for (auto elem_it = row.begin(); elem_it != row.end(); ++elem_it, ++V_it)
        {
            *elem_it = *V_it;
        }
    }
}

/**
 * @brief Fill a square matrix with the components of a vector
 *
 * @tparam ELEM_TYPE is the type of the vector components
 * @param V is the input vector
 * @return matrix<ELEM_TYPE>
 */
template <typename ELEM_TYPE>
matrix<ELEM_TYPE> gamma(const std::vector<ELEM_TYPE> &V)
{
    const size_t n = std::floor(std::sqrt(static_cast<double>(V.size())));

    if (n * n < V.size())
    {
        throw std::domain_error("The parameter cannot be reshaped in a square matrix");
    }

    matrix<ELEM_TYPE> M(n, std::vector<ELEM_TYPE>(n));

    fill_matrix(M, gamma);

    return M;
}

template<typename INTEGER_TYPE>
INTEGER_TYPE get_seed()
{
    using namespace Rcpp;

    Function sample_int("sample.int");

    Environment env = Environment::base_env();
    List machine = as<List>(env[".Machine"]);

    return as<INTEGER_TYPE>(sample_int(machine["integer.max"], 1));
}

template <typename SIGMA_VALUE_TYPE = double>
unsigned int significativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c,
                                    const size_t &n, const unsigned int m,
                                    const unsigned int &number_of_samples)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) <= c)
        {
            return 1;
        }

        return 0;
    };

    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    const size_t k = n * n;

    Rcpp::NumericVector V(k);

    const auto class_size = binom<size_t, mpz_class>(m + k - 1, m);

    unsigned int counter = 0;
    for (size_t i = 0; i < number_of_samples; ++i)
    {
        mpz_class rand = rand_generator.get_z_range(class_size);
        fill_iota(V, m, k, rand);
        Rcpp::NumericMatrix M(n, n, V.begin());
        counter += indicatorFunction(M);
    }

    return counter;
}


template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double significativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c, 
                              const size_t &n, const unsigned int m,
                              const unsigned int &number_of_samples)
{
    return static_cast<double>(significativityCounter(sigma, c, n, m, number_of_samples)) / number_of_samples;
}


template <typename SIGMA_VALUE_TYPE = double>
mpz_class significativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c,
                                 const size_t &n, const unsigned int m)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) <= c)
        {
            return 1;
        }

        return 0;
    };

    const size_t k = n * n;

    Rcpp::NumericVector V(k);

    const auto class_size = binom<size_t, mpz_class>(m + k - 1, m);

    mpz_class counter = 0;
    for (mpz_class i = 0; i < class_size; i = i+1)
    {
        fill_iota(V, m, k, i);
        Rcpp::NumericMatrix M(n, n, V.begin());
        counter = counter+indicatorFunction(M);
    }

    return counter;
}


template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double significativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c,
                              const size_t &n, const unsigned int m)
{
    mpq_class r(significativityCounter(sigma, c, n, m), binom<size_t, mpz_class>(m + n*n - 1, m));

    return r.get_d();
}


template <typename VECTOR_TYPE>
void fill_sample_probability_simplex(VECTOR_TYPE& V, gmp_randclass& rand_generator)
{
    double total{0.0};
    for (auto V_it = V.begin(); V_it != V.end(); ++V_it) {
        const mpf_class rand_value = rand_generator.get_f();
        *V_it = -std::log(rand_value.get_d());
        total += *V_it;
    }

    for (auto V_it = V.begin(); V_it != V.end(); ++V_it) {
        *V_it = *V_it / total;
    }
}


template <typename VECTOR_TYPE>
VECTOR_TYPE sample_probability_simplex(const size_t& k, gmp_randclass& rand_generator)
{
    VECTOR_TYPE V(k);

    fill_sample_probability_simplex(V, rand_generator);

    return V;
}


template <typename SIGMA_VALUE_TYPE = double>
unsigned int PsignificativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c, const size_t &n,
                                     const unsigned int &number_of_samples)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) <= c)
        {
            return 1;
        }

        return 0;
    };

    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    const size_t k = n * n;

    Rcpp::NumericVector V(k);

    unsigned int counter = 0;
    for (size_t i = 0; i < number_of_samples; ++i)
    {
        fill_sample_probability_simplex(V, rand_generator);
        Rcpp::NumericMatrix M(n, n, V.begin());
        counter += indicatorFunction(M);
    }

    return counter;
}

template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double Psignificativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c, const size_t &n,
                               const unsigned int &number_of_samples)
{
    return static_cast<double>(PsignificativityCounter(sigma, c, n, number_of_samples)) / number_of_samples;
}
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

    INDEX_TYPE C = binom<unsigned int, INDEX_TYPE>(m + k - 1, m);
    while (k > 1)
    {
        C = (k - 1) * C / (m + k - 1); // eval binom(m+(k-1)-2,m)

        unsigned int ml = m;
        while (m>0 && i>=C)
        {
            i = i - C;
            C = (m * C) / (m + k - 2); // eval binom((m-1)+k-2,(m-1))

            --m;
        }

        *ell_it = ml - m;
        ++ell_it;

        --k;
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
 * @return VECTOR_TYPE
 */
template <typename VECTOR_TYPE, typename INDEX_TYPE>
VECTOR_TYPE iota(const unsigned int m, const unsigned int k, const INDEX_TYPE i)
{
    VECTOR_TYPE ell(k);

    fill_iota(ell, m, k, i);

    return ell;
}

/**
 * @brief Get the random seed from the R environment
 * 
 * @tparam INTEGER_TYPE the random seed type
 * @return INTEGER_TYPE a integer 
 */
template<typename INTEGER_TYPE>
INTEGER_TYPE get_seed()
{
    using namespace Rcpp;

    Function sample_int("sample.int");

    Environment env = Environment::base_env();
    List machine = as<List>(env[".Machine"]);

    return as<INTEGER_TYPE>(sample_int(machine["integer.max"], 1));
}

/**
 * @brief Sample N times M_{n,m} and count how many samples M have sigma(M)<c
 * 
 * This function uniformly samples N times the set of the (n x n)-confusion
 * matrices whose elements sum up to m, i.e., M_{n,m}, and counts how many of 
 * the samples M are such that sigma(M)<c.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered confusion matrices
 * @param m is the number of tests used to build the considered confusion matrices
 * @param N is the number of confusion matrices to be sampled
 * @return the number of sampled confusion matrices M such that sigma(M)<c
 */
template <typename SIGMA_VALUE_TYPE = double>
unsigned int significativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c,
                                    const size_t &n, const unsigned int m,
                                    const unsigned int &N)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) < c)
        {
            return 1;
        }

        return 0;
    };

    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    const size_t k = n * n;

    Rcpp::NumericMatrix M(n);

    const auto class_size = binom<size_t, mpz_class>(m + k - 1, m);

    unsigned int counter = 0;
    for (size_t i = 0; i < N; ++i)
    {
        mpz_class rand = rand_generator.get_z_range(class_size);
        fill_iota(M, m, k, rand);
        counter += indicatorFunction(M);
    }

    return counter;
}

/**
 * @brief Sample N times M_{(s,.)} and count how many samples M have sigma(M)<c
 * 
 * This function uniformly samples N times the set of the (n x n)-confusion
 * matrices respecting s on rows, i.e., M_{(s,.)}, and counts how many of 
 * the samples M are such that sigma(M)<c.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param s is the vector of the classifications of the row-classifier
 * @param N is the number of confusion matrices to be sampled
 * @return the number of sampled confusion matrices M such that sigma(M)<c
 */
template <typename SIGMA_VALUE_TYPE = double>
unsigned int significativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c,
                                    const std::vector<unsigned int>& s,
                                    const unsigned int &N)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) < c)
        {
            return 1;
        }

        return 0;
    };

    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    const size_t n = s.size(); 
    Rcpp::NumericMatrix M(n);

    unsigned int counter = 0;
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < s.size(); ++j) {
            const auto& m = s[j];
            const auto class_size = binom<size_t, mpz_class>(m + n - 1, m);
            mpz_class rand = rand_generator.get_z_range(class_size);
            M(j, Rcpp::_) = iota<Rcpp::NumericVector>(m, n, rand);
        }
        counter += indicatorFunction(M);
    }

    return counter;
}

/**
 * @brief Estimate the sigma-significativity of c in M_{n,m}
 * 
 * This function estimates the sigma-significativity of c in M_{n,m} by using 
 * the Monte Carlo method.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered confusion matrices
 * @param m is the number of tests used to build the considered confusion matrices
 * @param N is the number of confusion matrices to be sampled
 * @return The Monte Carlo estimation of the sigma-significativity of c in M_{n,m}
 */
template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double significativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c,
                              const size_t &n, const unsigned int m,
                              const unsigned int &N)
{
    return static_cast<double>(significativityCounter(sigma, c, n, m, N)) / N;
}

/**
 * @brief Count how many of the matrices in M_{n,m} have sigma(M)<c
 * 
 * This function counts how many of the (n x n)-confusion matrices whose elements 
 * sum up to m, i.e., how many of the matrices M in M_{n,m}, are such that 
 * sigma(M)<c.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered confusion matrices
 * @param m is the number of tests used to build the considered confusion matrices
 * @return the number of confusion matrices M in M_{n,m} such that sigma(M)<c
 */
template <typename SIGMA_VALUE_TYPE = double>
mpz_class significativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c,
                                 const size_t &n, const unsigned int m)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) < c)
        {
            return 1;
        }

        return 0;
    };

    const size_t k = n * n;
    Rcpp::NumericMatrix M(n);

    const auto class_size = binom<size_t, mpz_class>(m + k - 1, m);

    mpz_class counter = 0;
    for (mpz_class i = 0; i < class_size; i = i+1)
    {
        fill_iota(M, m, k, i);
        counter = counter+indicatorFunction(M);
    }

    return counter;
}

/**
 * @brief Compute the sigma-significativity of c in M_{n,m}
 * 
 * This function computes the sigma-significativity of c in M_{n,m}.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered confusion matrices
 * @param m is the number of tests used to build the considered confusion matrices
 * @return The sigma-significativity of c in M_{n,m}
 */
template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double significativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c,
                              const size_t &n, const unsigned int m)
{
    mpq_class r(significativityCounter(sigma, c, n, m), binom<size_t, mpz_class>(m + n*n - 1, m));

    return r.get_d();
}

/**
 * @brief Sample the probability simplex
 * 
 * This function uniformly samples the probability simplex whose dimension is
 * `V.size()` and copies the sample in `V`.
 * 
 * @tparam VECTOR_TYPE is the sample output type
 * @param V is the vector in which the sample will be placed
 * @param rand_generator is a random generator
 */
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

/**
 * @brief Sample the k-dimensional probability simplex
 * 
 * This function uniformly samples the k-probability simplex.
 * 
 * @tparam VECTOR_TYPE is the sample output type
 * @param k is the dimension of the probability simplex to be sampled
 * @param rand_generator is a random generator
 */
template <typename VECTOR_TYPE>
VECTOR_TYPE sample_probability_simplex(const size_t& k, gmp_randclass& rand_generator)
{
    VECTOR_TYPE V(k);

    fill_sample_probability_simplex(V, rand_generator);

    return V;
}

/**
 * @brief Sample N times P_{n} and count how many samples M have sigma(M)<c
 * 
 * This function uniformly samples N times the set of the (n x n)-probability
 * matrices, i.e., P_{n}, and counts how many of the samples M are such that 
 * sigma(M)<c.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered probability matrices
 * @param N is the number of probability matrices to be sampled
 * @return the number of sampled probability matrices M such that sigma(M)<c
 */
template <typename SIGMA_VALUE_TYPE = double>
unsigned int PsignificativityCounter(const Rcpp::Function &sigma, const SIGMA_VALUE_TYPE &c, const size_t &n,
                                     const unsigned int &N)
{
    auto indicatorFunction = [&sigma, &c](const Rcpp::NumericMatrix &M)
    {
        if (Rcpp::as<double>(sigma(M)) < c)
        {
            return 1;
        }

        return 0;
    };

    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    Rcpp::NumericMatrix M(n);

    unsigned int counter = 0;
    for (size_t i = 0; i < N; ++i)
    {
        fill_sample_probability_simplex(M, rand_generator);
        counter += indicatorFunction(M);
    }

    return counter;
}

/**
 * @brief Estimate the sigma-significativity of c in P_{n}
 * 
 * This function estimates the sigma-significativity of c in P_{n} by using 
 * the Monte Carlo method.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param n is the number of rows/columns of the considered probability matrices
 * @param N is the number of probability matrices to be sampled
 * @return The Monte Carlo estimation of the sigma-significativity of c in P_{n}
 */
template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double Psignificativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c, const size_t &n,
                               const unsigned int &N)
{
    return static_cast<double>(PsignificativityCounter(sigma, c, n, N)) / N;
}

/**
 * @brief Estimate the sigma-significativity of c in M_{(s,.)}
 * 
 * This function estimates the sigma-significativity of c in M_{(s,.)} by using 
 * the Monte Carlo method.
 * 
 * @tparam SIGMA_VALUE_TYPE is the sigma-value type
 * @param sigma is an agreement measure
 * @param c is a sigma-value
 * @param s is the vector of the classifications of the row-classifier
 * @param N is the number of confusion matrices to be sampled
 * @return The Monte Carlo estimation of the sigma-significativity of c in M_{(s,.)}
 */
template <typename SIGMA_TYPE, typename SIGMA_VALUE_TYPE = double>
inline double significativity(const SIGMA_TYPE &sigma, const SIGMA_VALUE_TYPE &c,
                              const std::vector<unsigned int>& s,
                              const unsigned int &N)
{
    return static_cast<double>(significativityCounter(sigma, c, s, N)) / N;
}
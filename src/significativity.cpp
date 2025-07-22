#include <iostream>

#include <Rcpp.h>

#include "significativity.hpp"

double get_real(const SEXP &sexp)
{
    switch (TYPEOF(sexp))
    {
    case REALSXP:
    {
        return Rcpp::as<double>(sexp);
    }
    default:
        throw std::domain_error("The parameter is not a real value.");
    }
}

double get_unsigned_int(const SEXP &sexp)
{
    switch (TYPEOF(sexp))
    {
    case REALSXP:
    case INTSXP:
    {
        auto c_sexp = Rcpp::as<double>(sexp);

        if (c_sexp >= 0 && c_sexp == std::floor(c_sexp))
        {
            return static_cast<unsigned int>(c_sexp);
        }
    }
    default:
        throw std::domain_error("The parameter is not an unsigned integer value.");
    }
}

unsigned int validate_unsigned_int(const SEXP &value, const std::string &name)
{
    try
    {
        return get_unsigned_int(value);
    }
    catch (std::domain_error &e)
    {
        ::Rf_error("The parameter \"%s\" must be an unsigned positive integer.", name.c_str());
    }
}

double validate_double(const SEXP &value, const std::string &name)
{
    try
    {
        return get_real(value);
    }
    catch (std::domain_error &e)
    {
        ::Rf_error("The parameter \"%s\" must be a real value.", name.c_str());
    }
}

double M_significativity(const Rcpp::Function &sigma, const SEXP &c, const SEXP &n,
                         const SEXP m, const SEXP &number_of_samples)
{
    auto c_c = validate_double(c, "c");
    auto c_n = validate_unsigned_int(n, "n");
    auto c_m = validate_unsigned_int(m, "m");

    if (number_of_samples == R_NilValue)
    {
        return significativity(sigma, c_c, c_n, c_m);
    }
    auto c_number_of_samples = validate_unsigned_int(number_of_samples, "number_of_samples");

    return significativity(sigma, c_c, c_n, c_m, c_number_of_samples);
}

double P_significativity(const Rcpp::Function &sigma, const SEXP &c, const SEXP &n,
                         const SEXP &number_of_samples)
{
    auto c_c = validate_double(c, "c");
    auto c_n = validate_unsigned_int(n, "n");
    auto c_number_of_samples = validate_unsigned_int(number_of_samples, "number_of_samples");

    return Psignificativity(sigma, c_c, c_n, c_number_of_samples);
}

double MP_significativity(const Rcpp::Function &sigma, const SEXP &c, const SEXP &n,
                          const SEXP m, const SEXP &number_of_samples)
{
    if (m == R_NilValue)
    {
        return P_significativity(sigma, c, n, number_of_samples);
    }
    return M_significativity(sigma, c, n, m, number_of_samples);
}

std::vector<unsigned int> validate_unsigned_int_vect(const SEXP &value, const std::string &name)
{
    if (! Rcpp::is<Rcpp::IntegerVector>(value)) {
        ::Rf_error("The parameter \"%s\" must be a vector of unsigned positive integers.", name.c_str());
    }

    return Rcpp::as<std::vector<unsigned int>>(value);
}

double marginal_significativity(const Rcpp::Function &sigma, const SEXP &c,
                                const SEXP &s, const SEXP &number_of_samples)
{
    auto c_c = validate_double(c, "c");
    auto c_s = validate_unsigned_int_vect(s, "s");
    auto c_number_of_samples = validate_unsigned_int(number_of_samples, "number_of_samples");

    return significativity(sigma, c_c, c_s, c_number_of_samples);
}

Rcpp::NumericVector sample_prob_simplex(const size_t& k)
{
    gmp_randclass rand_generator(gmp_randinit_default);
    rand_generator.seed(get_seed<int>());

    return sample_probability_simplex<Rcpp::NumericVector>(k, rand_generator);
}

Rcpp::NumericVector weak_composition_enum(const unsigned int& m, const unsigned int& k, const size_t index)
{
    return iota<Rcpp::NumericVector>(m, k, index);
}

using namespace Rcpp;

RCPP_MODULE(rSignificativity)
{

//' @name significativity
//' @title The \eqn{\sigma}-significativity of \eqn{c}
//' @description This function evaluates the \eqn{\sigma}-significativity of an
//'   agreement value \eqn{c} in either \eqn{\mathcal{M}_{n,m}} or
//'   \eqn{\mathcal{P}_{n}} depending on the call parameters.
//'
//'   When `m` is a natural number, the \eqn{\sigma}-significativity of \eqn{c}
//'   in \eqn{\mathcal{M}_{n,m}} is evaluated. If `number_of_samples` is
//'   a natural number, the significativity is estimated by using the Monte
//'   Carlo method. Otherwise, the computation considers all the
//'   \eqn{n \times n}-confusion matrices.
//'
//'   When instead `m` is set to `NULL`, the function uses the Monte Carlo
//'   method to estimate the \eqn{\sigma}-significativity of \eqn{c} in
//'   \eqn{\mathcal{P}_{n}}.
//' @param sigma An agreement measure.
//' @param c An agreement value.
//' @param n The number of rows/columns of the confusion matrix.
//' @param m The sum of the confusion matrix elements. When set to `NULL`,
//'   the function estimate the \eqn{\sigma}-significativity of \eqn{c} in
//'   \eqn{\mathcal{P}_{n}} (default: `NULL`).
//' @param number_of_samples The number of samples used to evaluate the
//'   significativity by using Monte Carlo method (default: 10000).
//' @return When `m` is a natural number, the function returns the
//'   \eqn{\sigma}-significativity of \eqn{c} in \eqn{\mathcal{M}_{n,m}}.
//'   If instead `m` is set to `NULL`, the \eqn{\sigma}-significativity
//' of \eqn{c} in \eqn{\mathcal{P}_{n}}.
//' @examples
//' # evaluate kappa-significativity of 0.5 in M_{2,5} with 10000 samples
//' significativity(cohen_kappa, 0.5, 2, 5)
//'
//' # evaluate kappa-significativity of 0.5 in M_{2,5} with 1000 samples
//' significativity(cohen_kappa, 0.5, 2, 5, number_of_samples = 1000)
//'
//' # exactly compute kappa-significativity of 0.5 in M_{2,5}
//' significativity(cohen_kappa, 0.5, 2, 5, number_of_samples = NULL)
//'
//' # evaluate kappa-significativity of 0.5 in P_{2} with 1000 samples
//' significativity(cohen_kappa, 0.5, 2, number_of_samples = 1000)
//'
//' # evaluate kappa-significativity of 0.5 in P_{2} with 10000 samples
//' significativity(cohen_kappa, 0.5, 2)
//'
//' # successive calls to Monte Carlo methods may produce different results
//' significativity(cohen_kappa, 0.5, 2)
//'
//' # setting the random seed before the call guarantee repeatability
//' set.seed(1)
//' significativity(cohen_kappa, 0.5, 2)
//'
//' set.seed(1)
//' significativity(cohen_kappa, 0.5, 2)
    function("significativity", &MP_significativity,
             List::create(_["sigma"], _["c"], _["n"], _["m"] = R_NilValue,
                          _["number_of_samples"] = 10000),
             "Estimate the sigma-significativity of c in P_{n}");

//' @name marginal.significativity
//' @title Estimate the marginal \eqn{\sigma}-significativity of \eqn{c}
//' @description This function evaluates the \eqn{\sigma}-significativity
//'   of an agreement value \eqn{c} in \eqn{\mathcal{M}_{(s,\cdot)}}.
//'
//'   \eqn{\mathcal{M}_{(s,\cdot)}} is the set of the confusion matrices
//'   having \eqn{|s|} rows and columns such that the elements in the 
//'   \eqn{i}-th row sum up to \eqn{s(i)}.
//'
//'   The \eqn{\sigma}-significativity of \eqn{c} in
//'   \eqn{\mathcal{M}_{(s,\cdot)}} is the ratio between the number of
//'   matrices \eqn{M \in \mathcal{M}_{(s,\cdot)}} such that
//'   \eqn{\sigma(M) < c} and the cardinality of
//'   \eqn{\mathcal{M}_{(s,\cdot)}}. This corresponds to the
//'   probability of choosing with uniform distribution a matrix
//'   \eqn{M \in \mathcal{M}_{(s,\cdot)}} such that \eqn{\sigma(M) < c}.
//' @param sigma An agreement measure.
//' @param c An agreement value.
//' @param s A vector of natural values.
//' @param number_of_samples The number of samples used to evaluate the
//'   significativity by using Monte Carlo method (default: 10000).
//' @return The function returns the \eqn{\sigma}-significativity of
//'   \eqn{c} in \eqn{\mathcal{M}_{(s,\cdot)}} as estimated by the
//'   Monte Carlo method using `number_of_samples` samples.
//' @examples
//' # define a vector defining the number of classifications of a
//' # two-class classifier 
//' s <- as.integer(c(10, 1))
//'
//' # gauge kappa-significativity of 0.5 in M_{(s,.)} with 10000 samples
//' marginal.significativity(cohen_kappa, 0.5, s)
//'
//' # define a different vector defining the number of classifications of a
//' # two-class classifier 
//' s2 <- as.integer(c(6, 5))
//'
//' # define a vector defining the number of classifications of a five-class
//' # classifier 
//' s3 <- as.integer(c(6, 5, 3, 8, 5))
//'
//' # gauge kappa-significativity of 0.5 in M_{(s3,.)} with 10000 samples
//' marginal.significativity(cohen_kappa, 0.5, s3)
//'
//' # gauge kappa-significativity of 0.5 in M_{(s3,.)} with 40000 samples
//' marginal.significativity(cohen_kappa, 0.5, s3, number_of_samples=40000)
//'
//' # gauge kappa-significativity of 0.5 in M_{(s2,.)} with 10000 samples
//' marginal.significativity(cohen_kappa, 0.5, s2)
//'
//' # define a (2x2)-confusion matrix
//' M <- matrix(as.integer(c(9, 0, 1, 1)), nrow=2)
//' M
//' cohen_kappa(M)
//'
//' # gauge kappa-significativity of cohen_kappa(M) in M_{(rowSums(M),.)}
//' # with 10000 samples
//' marginal.significativity(cohen_kappa, cohen_kappa(M),
//'                          as.integer(rowSums(M)))
//'
//' # define another (2x2)-confusion matrix
//' M2 <- matrix(as.integer(c(5, 0, 1, 5)), nrow=2)
//' M2
//' cohen_kappa(M2)
//'
//' # gauge kappa-significativity of cohen_kappa(M2) in M_{(rowSums(M2),.)}
//' # with 10000 samples
//' marginal.significativity(cohen_kappa, cohen_kappa(M2),
//'                          as.integer(rowSums(M2)))
    function("marginal.significativity", &marginal_significativity,
             List::create(_["sigma"], _["c"], _["s"],
                          _["number_of_samples"] = 10000),
             "Estimate the sigma-significativity of c in M_{(s,.)}");

//' @name weak_composition_enum
//' @title The lexicographic weak composition enumeration
//' @description  This function implements the lexicographic order enumerator 
//'   for weak compositions of \eqn{m} in \eqn{k} parts.
//'
//'   A weak composition of \eqn{m} in \eqn{k} parts is a vector of \eqn{k}
//'   natural numbers whose sum is \eqn{m}. \eqn{\mathcal{C}_{m,k}} is the
//'   set of all the weak composition of \eqn{m} in \eqn{k} parts.
//'
//'   The lexicographic order \eqn{<_l} over \eqn{\mathcal{C}_{m,k}} is an
//'   order between weak compositions of \eqn{m} in \eqn{k} parts such that
//'   \eqn{\langle a_1,\ldots a_k \rangle <_l \langle b_1,\ldots b_k \rangle}
//'   when there exist \eqn{i \in [1, k]} such that for all
//'   \eqn{j \in [1, i-1]} it holds that \eqn{a_j = b_j} and \eqn{a_i<b_i}.
//' @param m The sum of the natural values in the output vector.
//' @param k The size of the output vector.
//' @param index The index of the output vector in the lexicographic order.
//' @return The `index`-th weak composition in the lexicographic order of
//'   \eqn{\mathcal{C}_{m,k}}.
//' @examples
//' # evaluate the 1st element in the lexicographic order of C_{100, 5}
//' weak_composition_enum(100, 5, 0)
//'
//' # evaluate it again
//' weak_composition_enum(100, 5, 0)
//'
//' # evaluate the 42137th element in the lexicographic order of C_{100, 5}
//' weak_composition_enum(100, 5, 42137)
    function("weak_composition_enum", &weak_composition_enum,
             List::create(_["m"], _["k"], _["index"]),
             "Enumerator function for the weak compositions in C_{m,k}.");

//' @name sample_prob_simplex
//' @title Uniformly sampling \eqn{\Delta^{(k-1)}}
//' @description This function samples the k-dimensional probability simplex
//'   \eqn{\Delta^{(k-1)}=\left\{ \langle x_1, \ldots x_k \rangle \in \mathbb{R}_{\geq 0}^k\, |\, \sum_{i=1}^k x_i = 1\right\}}
//'   using the uniform probability distribution.
//' @param k The dimension of the probability simplex to be sampled.
//' @return A sample of \eqn{\Delta^{(k-1)}} taken according the uniform
//'   probability distribution.
//' @examples
//' # sample the 5-dimensional probability simplex
//' sample_prob_simplex(5)
//'
//' # successive calls may produce different results
//' sample_prob_simplex(5)
//'
//' # setting the random seed before the call guarantee repeatability
//' set.seed(1)
//' sample_prob_simplex(5)
//'
//' set.seed(1)
//' sample_prob_simplex(5)
    function("sample_prob_simplex", &sample_prob_simplex,
             List::create(_["k"]),
             "Sample the k-dimensional probability simplex");
}
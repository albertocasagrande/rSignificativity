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

using namespace Rcpp;

RCPP_MODULE(rSignificativity)
{

//' @name significativity
//' @title The \eqn{\sigma}-significativity of \eqn{c}
//' @description This function evaluates the significativity of a
//'   \eqn{\sigma}-value in either \eqn{\mathcal{M}_{n,m}} or
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
//' @param sigma A statistical coefficient.
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
}
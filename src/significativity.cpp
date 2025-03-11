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

double M_significativity(const Rcpp::Function &sigma, const SEXP &n, const SEXP m,
                         const SEXP &c, const SEXP &number_of_samples)
{
    auto c_n = validate_unsigned_int(n, "n");
    auto c_m = validate_unsigned_int(m, "m");
    auto c_c = validate_double(c, "c");

    if (number_of_samples == R_NilValue)
    {
        return significativity(sigma, c_n, c_m, c_c);
    }
    auto c_number_of_samples = validate_unsigned_int(number_of_samples, "number_of_samples");

    return significativity(sigma, c_n, c_m, c_c, c_number_of_samples);
}

double P_significativity(const Rcpp::Function &sigma, const SEXP &n,
                         const SEXP &c, const SEXP &number_of_samples)
{
    auto c_n = validate_unsigned_int(n, "n");
    auto c_c = validate_double(c, "c");
    auto c_number_of_samples = validate_unsigned_int(number_of_samples, "number_of_samples");

    return Psignificativity(sigma, c_n, c_c, c_number_of_samples);
}

using namespace Rcpp;

RCPP_MODULE(rSignificativity)
{
    //' @name M_significativity
    //' @title Estimate the $\sigma$-significativity of $c$ in $\mathcal{M}_{n,m}$
    //' @description This method evaluates the significativity of a sigma-value. It
    //'   uses the Monte Carlo method by default and uniformly samples the set of all
    //'   the $n \times n$-confusion matrices summing up to $m$. When
    //'   `number_of_samples` is set to `NULL`, the computation considers all the
    //'   $n \times n$-confusion matrices.
    //' @param sigma A statistical coefficient.
    //' @param n The number of rows/columns of the confusion matrix.
    //' @param m The sum of the confusion matrix elements.
    //' @param number_of_samples The number of samples used to evaluate the
    //'   significativity by using Monte Carlo method (default: 10000).
    //'
    //'   When `number_of_samples` is set to `NULL`, the function don't use the Monte
    //'   Carlo method and computes the significativity by considering all the
    //'   confusion matrices in $\mathcal{M}_{n,m}$
    //' @return The $\sigma$-significativity of $c$ in $\mathcal{M}_{n,m}$.
    function("M_significativity", &M_significativity,
             List::create(_["sigma"], _["n"], _["m"], _["c"], _["number_of_samples"] = 10000),
             "Estimate the sigma-significativity of c in M_{n,m}");

    //' @name P_significativity
    //' @title Estimate the $\sigma$-significativity of $c$ in $\mathcal{P}_{n}$
    //' @description This method evaluates the significativity of a sigma-value. It
    //'   uses the Monte Carlo method and uniformly samples the set of all
    //'   the $n \times n$-probability matrices.
    //' @param sigma A statistical coefficient.
    //' @param n The number of rows/columns of the confusion matrix.
    //' @param number_of_samples The number of samples used to evaluate the
    //'   significativity by using Monte Carlo method (default: 10000).
    //' @return The $\sigma$-significativity of $c$ in $\mathcal{P}_{n}$.
    function("P_significativity", &P_significativity,
             List::create(_["sigma"], _["n"], _["c"], _["number_of_samples"] = 10000),
             "Estimate the sigma-significativity of c in P_{n}");
}
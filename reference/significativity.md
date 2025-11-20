# The \\\sigma\\-significativity of \\c\\

This function evaluates the \\\sigma\\-significativity of an agreement
value \\c\\ in either \\\mathcal{M}\_{n,m}\\ or \\\mathcal{P}\_{n}\\
depending on the call parameters.

When `m` is a natural number, the \\\sigma\\-significativity of \\c\\ in
\\\mathcal{M}\_{n,m}\\ is evaluated. If `number_of_samples` is a natural
number, the significativity is estimated by using the Monte Carlo
method. Otherwise, the computation considers all the \\n \times
n\\-confusion matrices.

When instead `m` is set to `NULL`, the function uses the Monte Carlo
method to estimate the \\\sigma\\-significativity of \\c\\ in
\\\mathcal{P}\_{n}\\.

## Arguments

- sigma:

  An agreement measure.

- c:

  An agreement value.

- n:

  The number of rows/columns of the confusion matrix.

- m:

  The sum of the confusion matrix elements. When set to `NULL`, the
  function estimate the \\\sigma\\-significativity of \\c\\ in
  \\\mathcal{P}\_{n}\\ (default: `NULL`).

- number_of_samples:

  The number of samples used to evaluate the significativity by using
  Monte Carlo method (default: 10000).

## Value

When `m` is a natural number, the function returns the
\\\sigma\\-significativity of \\c\\ in \\\mathcal{M}\_{n,m}\\. If
instead `m` is set to `NULL`, the \\\sigma\\-significativity of \\c\\ in
\\\mathcal{P}\_{n}\\.

## Examples

``` r
# evaluate kappa-significativity of 0.5 in M_{2,5} with 10000 samples
significativity(cohen_kappa, 0.5, 2, 5)
#> [1] 0.8263

# evaluate kappa-significativity of 0.5 in M_{2,5} with 1000 samples
significativity(cohen_kappa, 0.5, 2, 5, number_of_samples = 1000)
#> [1] 0.804

# exactly compute kappa-significativity of 0.5 in M_{2,5}
significativity(cohen_kappa, 0.5, 2, 5, number_of_samples = NULL)
#> [1] 0.8214286

# evaluate kappa-significativity of 0.5 in P_{2} with 1000 samples
significativity(cohen_kappa, 0.5, 2, number_of_samples = 1000)
#> [1] 0.886

# evaluate kappa-significativity of 0.5 in P_{2} with 10000 samples
significativity(cohen_kappa, 0.5, 2)
#> [1] 0.897

# successive calls to Monte Carlo methods may produce different results
significativity(cohen_kappa, 0.5, 2)
#> [1] 0.8982

# setting the random seed before the call guarantee repeatability
set.seed(1)
significativity(cohen_kappa, 0.5, 2)
#> [1] 0.9012

set.seed(1)
significativity(cohen_kappa, 0.5, 2)
#> [1] 0.9012
```

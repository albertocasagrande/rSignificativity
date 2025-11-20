# Estimate the marginal \\\sigma\\-significativity of \\c\\

This function evaluates the \\\sigma\\-significativity of an agreement
value \\c\\ in \\\mathcal{M}\_{(s,\cdot)}\\.

\\\mathcal{M}\_{(s,\cdot)}\\ is the set of the confusion matrices having
\\\|s\|\\ rows and columns such that the elements in the \\i\\-th row
sum up to \\s(i)\\.

The \\\sigma\\-significativity of \\c\\ in \\\mathcal{M}\_{(s,\cdot)}\\
is the ratio between the number of matrices \\M \in
\mathcal{M}\_{(s,\cdot)}\\ such that \\\sigma(M) \< c\\ and the
cardinality of \\\mathcal{M}\_{(s,\cdot)}\\. This corresponds to the
probability of choosing with uniform distribution a matrix \\M \in
\mathcal{M}\_{(s,\cdot)}\\ such that \\\sigma(M) \< c\\.

## Arguments

- sigma:

  An agreement measure.

- c:

  An agreement value.

- s:

  A vector of natural values.

- number_of_samples:

  The number of samples used to evaluate the significativity by using
  Monte Carlo method (default: 10000).

## Value

The function returns the \\\sigma\\-significativity of \\c\\ in
\\\mathcal{M}\_{(s,\cdot)}\\ as estimated by the Monte Carlo method
using `number_of_samples` samples.

## Examples

``` r
# define a vector defining the number of classifications of a
# two-class classifier 
s <- as.integer(c(10, 1))

# gauge kappa-significativity of 0.5 in M_{(s,.)} with 10000 samples
marginal.significativity(cohen_kappa, 0.5, s)
#> [1] 0.9032

# define a different vector defining the number of classifications of a
# two-class classifier 
s2 <- as.integer(c(6, 5))

# define a vector defining the number of classifications of a five-class
# classifier 
s3 <- as.integer(c(6, 5, 3, 8, 5))

# gauge kappa-significativity of 0.5 in M_{(s3,.)} with 10000 samples
marginal.significativity(cohen_kappa, 0.5, s3)
#> [1] 0.9995

# gauge kappa-significativity of 0.5 in M_{(s3,.)} with 40000 samples
marginal.significativity(cohen_kappa, 0.5, s3, number_of_samples=40000)
#> [1] 0.99955

# gauge kappa-significativity of 0.5 in M_{(s2,.)} with 10000 samples
marginal.significativity(cohen_kappa, 0.5, s2)
#> [1] 0.8552

# define a (2x2)-confusion matrix
M <- matrix(as.integer(c(9, 0, 1, 1)), nrow=2)
M
#>      [,1] [,2]
#> [1,]    9    1
#> [2,]    0    1
cohen_kappa(M)
#> [1] 0.6206897

# gauge kappa-significativity of cohen_kappa(M) in M_{(rowSums(M),.)}
# with 10000 samples
marginal.significativity(cohen_kappa, cohen_kappa(M),
                         as.integer(rowSums(M)))
#> [1] 0.914

# define another (2x2)-confusion matrix
M2 <- matrix(as.integer(c(5, 0, 1, 5)), nrow=2)
M2
#>      [,1] [,2]
#> [1,]    5    1
#> [2,]    0    5
cohen_kappa(M2)
#> [1] 0.8196721

# gauge kappa-significativity of cohen_kappa(M2) in M_{(rowSums(M2),.)}
# with 10000 samples
marginal.significativity(cohen_kappa, cohen_kappa(M2),
                         as.integer(rowSums(M2)))
#> [1] 0.9548
```

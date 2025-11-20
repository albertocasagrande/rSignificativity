# Uniformly sampling \\\Delta^{(k-1)}\\

This function samples the k-dimensional probability simplex
\\\Delta^{(k-1)}=\left\\ \langle x_1, \ldots x_k \rangle \in
\mathbb{R}\_{\geq 0}^k\\ \|\\ \sum\_{i=1}^k x_i = 1\right\\\\ using the
uniform probability distribution.

## Arguments

- k:

  The dimension of the probability simplex to be sampled.

## Value

A sample of \\\Delta^{(k-1)}\\ taken according the uniform probability
distribution.

## Examples

``` r
# sample the 5-dimensional probability simplex
sample_prob_simplex(5)
#> [1] 0.1701992 0.1082816 0.2597004 0.2710165 0.1908022

# successive calls may produce different results
sample_prob_simplex(5)
#> [1] 0.01736597 0.16952667 0.43687385 0.18650819 0.18972533

# setting the random seed before the call guarantee repeatability
set.seed(1)
sample_prob_simplex(5)
#> [1] 0.05830519 0.15783181 0.55531351 0.01952378 0.20902571

set.seed(1)
sample_prob_simplex(5)
#> [1] 0.05830519 0.15783181 0.55531351 0.01952378 0.20902571
```

# Agreement Value Significativity

Agreement measures, such as Cohen’s $`\kappa`$(Cohen 1960) and
$`\textrm{IA}`$(Casagrande et al. 2020), gauge the agreement between two
discrete classifiers mapping their confusion matrices into real values.

``` r

# create a confusion matrix
M <- matrix(1:9, nrow = 3, ncol = 3)
print(M)
#>      [,1] [,2] [,3]
#> [1,]    1    4    7
#> [2,]    2    5    8
#> [3,]    3    6    9

# rSignificativity implements some agreement measures
rSignificativity::cohen_kappa(M)
#> [1] -0.04166667
rSignificativity::scott_pi(M)
#> [1] -0.05633803
rSignificativity::IA(M)
#> [1] 0.005631984
```

When a golden standard classifier is available, the best among a set of
classifiers can be selected by picking the one whose confusion matrix
with respect to the golden standard has the maximal agreement value.

The agreement values alone lack of a meaningfulness indication.
Understanding whether two classifiers *significantly* agree on the
considered dataset is not possible by simply looking at the
corresponding agreement value.

Let $`\mathcal{M}_{n,m}`$ be the set of $`n \times n`$-confusion
matrices built over a dataset having size $`m`$, i.e., the elements of
any matrix $`\mathcal{M}_{n,m}`$ sum up to $`m`$. For any agreement
measure $`\sigma`$, the $`\sigma`$-significativity of $`c`$ in
$`\mathcal{M}_{n,m}`$(Casagrande et al. 2025) is \$\$ \varrho\_{\sigma,
n, m}(c) \stackrel{\tiny\textrm{def}}{=} \frac{\left\|\left\\M \in
\mathcal{M}\_{n,m} \\ \|\\ \sigma(M) \<
c\right\\\right\|}{\left\|\mathcal{M}\_{n,m}\right\|}. \$\$ The
$`\sigma`$-significativity is the \[0,1\]-normalised number of matrices
in $`\mathcal{M}_{n,m}`$ having a $`\sigma`$-value greater than $`c`$.
From a different point of view, it corresponds to the probability of
choosing by chance a matrix whose $`\sigma`$-value greater than $`c`$.

Since $`|\mathcal{M}_{n,m}|=\binom{n^2+m-1}{m}`$, computing
$`\varrho_{\sigma, n, m}(c)`$ takes time
$`O\left(n^2\binom{n^2+m-1}{m}\right)`$. Nevertheless, the Monte Carlo
method (Metropolis and Ulam 1949) can estimate
$`\varrho_{\sigma, n, m}(c)`$ with $`N`$ samples in time
$`\Theta(N(n^2+m))`$ with an error proportional to $`1/\sqrt{N}`$ (e.g.,
see (Mackay 1998)).

The function
[`significativity()`](https://albertocasagrande.github.io/rSignificativity/reference/significativity.md)
lets us to both exactly compute $`\varrho_{\sigma, n, m}(c)`$ and
estimate it by using the Monte Carlo method.

``` r

library("rSignificativity")

# evaluate kappa-significativity of 0.5 in M_{2,100} with 10000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 100)
#> [1] 0.888

# evaluate kappa-significativity of 0.5 in M_{2,100} with 1000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 100, N = 1000)
#> [1] 0.897

# evaluate kappa-significativity of 0.5 in M_{2,5} with 1000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 5, N = 1000)
#> [1] 0.835

# exactly compute kappa-significativity of 0.5 in M_{2,5}
significativity(cohen_kappa, c = 0.5, n = 2, m = 5, N = NULL)
#> [1] 0.8214286
```

The $`\sigma`$-significativity of $`c`$ in
$`\mathcal{P}_{n}`$(Casagrande et al. 2025), where $`\mathcal{P}_{n}`$
is the set of all the $`n \times n`$-probability matrices, is \$\$
\rho\_{\sigma, n}(c) \stackrel{\tiny\textrm{def}}{=}
\frac{V\left(\left\\M \in \mathcal{P}\_{n} \\ \|\\ \sigma(M) \<
c\right\\\right)}{V(\Delta^{(n^2-1)})} \$\$ where
$`\Delta^{(k-1)} =\left\{ \langle x_1, \ldots x_k \rangle \in \mathbb{R}_{\geq 0}^k\, |\, \sum_{i=1}^k x_i = 1\right\}`$
is the $`k`$-dimensional probability simplex and $`V(\cdot)`$ is the
$`n^2-1`$-dimensional Lebesgue measure (Casagrande et al. 2025).

If $`\sigma(M) < c`$ is definable in an o-minimal theory (van den Dries
1998), such as in the case of Cohen’s $`\kappa`$ and $`\textrm{IA}*`$,
then
``` math
\lim_{m \rightarrow +\infty}\varrho_{\sigma, n, m}(c) = \rho_{\sigma, n}(c).
```

We can estimate $`\rho_{\sigma, n}(c)`$ by using the Monte Carlo method
with $`N`$ samples in time $`\Theta(n^2 N)`$. The function
[`significativity()`](https://albertocasagrande.github.io/rSignificativity/reference/significativity.md)
also implements this algorithm.

``` r

# evaluate kappa-significativity of 0.5 in P_{2} with 1000 samples
significativity(cohen_kappa, 0.5, n = 2, N = 1000)
#> [1] 0.895

# evaluate kappa-significativity of 0.5 in P_{2} with 10000 samples
significativity(cohen_kappa, 0.5, n = 2)
#> [1] 0.8979

# successive calls to Monte Carlo methods may produce different results
significativity(cohen_kappa, 0.5, n = 2)
#> [1] 0.8961

# setting the random seed before the call guarantee repeatability
set.seed(1)
significativity(cohen_kappa, 0.5, n = 2)
#> [1] 0.9012

set.seed(1)
significativity(cohen_kappa, 0.5, n = 2)
#> [1] 0.9012
```

## Ground-truth and Agreement: Marginal Significativity

The $`\sigma`$-significativity does not distinguish between the two
classifiers that produce a confusion matrix: both of them are assumed to
be equally relevant and any matrix in $`\mathcal{M}_{n,m}`$ is
admissible. However, in many contexts one of the two classifiers is the
gold standard, or the *de facto* standard, of the domain. Its
classifications are correct by definition and they should not change:
only the other classifier may vary.

For instance, if the gold standard labels 3 entries as $`A`$ and 17
entries as $`B`$ in a dataset of size 20, the only admissible confusion
matrices are those whose rows associated with $`A`$ and $`B`$ sum up to
3 and 17, respectively. Not all the matrices in $`\mathcal{M}_{n,m}`$
are admissible anymore and the significativity of an agreement value
should take this into account. The *marginal $`\sigma`$-significativity*
(Casagrande and Pagliarini 2025) addresses this issue.

### Definition

An $`n \times n`$ confusion matrix $`M`$*respects*
$`\mathbf{s} \in \mathbb{N}^n`$ on rows when $`\mathbf{s}(i)`$ is the
sum of the elements in the $`i`$-th row of $`M`$ for all
$`i \in [1, n]`$. The set of the confusion matrices respecting
$`\mathbf{s}`$ on rows is \$\$ \mathcal{M}\_{(\mathbf{s},\cdot)}
\stackrel{\tiny\textrm{def}}{=} \left\\M \in \mathbb{N}^{n\times n}
\\\Big\|\\ \mathbf{s}(j) = \sum\_{i=1}^n M(j,i) \textrm{ for all } j \in
\[1,n\]\right\\. \$\$ The set $`\mathcal{M}_{(\cdot,\mathbf{s})}`$ of
the confusion matrices respecting $`\mathbf{s}`$ on columns is defined
analogously. Whenever $`m = \sum_{i=1}^n \mathbf{s}(i)`$, both
$`\mathcal{M}_{(\mathbf{s},\cdot)}`$ and
$`\mathcal{M}_{(\cdot,\mathbf{s})}`$ are subsets of
$`\mathcal{M}_{n,m}`$.

The *marginal $`\sigma`$-significativity of $`c`$ respecting
$`\mathbf{s}`$* is \$\$ \varrho\_{\sigma, (\mathbf{s},\cdot)}(c)
\stackrel{\tiny\textrm{def}}{=} \frac{\left\|\left\\M \in
\mathcal{M}\_{(\mathbf{s},\cdot)} \\ \|\\ \sigma(M) \<
c\right\\\right\|}{\left\|\mathcal{M}\_{(\mathbf{s},\cdot)}\right\|}.
\$\$ It is the probability of selecting by chance, with uniform
distribution, a confusion matrix whose rows sum up to $`\mathbf{s}`$ and
whose $`\sigma`$-value is lower than $`c`$. If the rows of the matrix
represent the gold standard, $`\varrho_{\sigma, (\mathbf{s},\cdot)}(c)`$
is the probability that a classifier randomly labelling the dataset
entries reaches an agreement value lower than $`c`$ with the gold
standard.

The $`j`$-th row of every matrix in $`\mathcal{M}_{(\mathbf{s},\cdot)}`$
is a *weak composition* of $`\mathbf{s}(j)`$ in $`n`$ parts, i.e., a
tuple of $`n`$ natural numbers summing up to $`\mathbf{s}(j)`$. The set
$`\mathcal{C}_{m,k}`$ of the weak compositions of $`m`$ in $`k`$ parts
has cardinality $`\binom{m+k-1}{m}`$(Bóna 2016). Since the rows of the
matrices in $`\mathcal{M}_{(\mathbf{s},\cdot)}`$ are pairwise
independent,
``` math
\left|\mathcal{M}_{(\mathbf{s},\cdot)}\right| = \prod_{j=1}^n \binom{\mathbf{s}(j)+n-1}{\mathbf{s}(j)}.
```
Hence, enumerating $`\mathcal{M}_{(\mathbf{s},\cdot)}`$ takes time
$`O\left(n^2\binom{n^2+m-1}{m}\right)`$. Nevertheless, the Monte Carlo
method (Metropolis and Ulam 1949) can estimate
$`\varrho_{\sigma, (\mathbf{s},\cdot)}(c)`$ using $`N`$ samples with an
error in $`O(1/\sqrt{N})`$ in time $`O(N(n^2+m)) \cap \Omega(N n^2)`$.

This is the algorithm implemented by
[`marginal.significativity()`](https://albertocasagrande.github.io/rSignificativity/reference/marginal.significativity.md),
whose parameter `N` sets the number of samples (default: 10000).

``` r

s <- as.integer(c(10, 20, 30))

set.seed(1)
marginal.significativity(cohen_kappa, 0.2, s, N = 1000)
#> [1] 0.799

set.seed(1)
marginal.significativity(cohen_kappa, 0.2, s, N = 100000)
#> [1] 0.8186
```

## References

Bóna, Miklós. 2016. *A Walk Through Combinatorics: An Introduction to
Enumeration and Graph Theory, Second Edition*. 4th ed. World Scientific
Publishing Company. <https://doi.org/10.1142/10258>.

Casagrande, Alberto, Francesco Fabris, and Rossano Girometti. 2020.
“Extending Information Agreement by Continuity.” *2020 IEEE
International Conference on Bioinformatics and Biomedicine (BIBM)*,
1432–39. <https://doi.org/10.1109/BIBM49941.2020.9313173>.

Casagrande, Alberto, Francesco Fabris, Rossano Girometti, and Roberto
Pagliarini. 2025. “Significativity Indices for Agreement Values.”
*Statistics and Computing* 35 (6): 197.
<https://doi.org/10.1007/s11222-025-10728-1>.

Casagrande, Alberto, and Roberto Pagliarini. 2025. “A Marginal
Significativity Index for Agreement Values.” *2025 IEEE International
Conference on Bioinformatics and Biomedicine (BIBM)*, 5392–99.
<https://doi.org/10.1109/BIBM66473.2025.11356859>.

Cohen, Jacob. 1960. “A Coefficient of Agreement for Nominal Scales.”
*Educational and Psychological Measurement* 20: 37–46.
<https://doi.org/10.1177/001316446002000104>.

Mackay, D. J. C. 1998. “Introduction to Monte Carlo Methods.” In
*Learning in Graphical Models*, edited by Michael I. Jordan. Springer
Netherlands. <https://doi.org/10.1007/978-94-011-5014-9_7>.

Metropolis, Nicholas, and Stanisław Ulam. 1949. “The Monte Carlo
Method.” *Journal of the American Statistical Association* 44 (247):
335–41. <https://doi.org/10.1080/01621459.1949.10483310>.

van den Dries, Lou P. D. 1998. *Tame Topology and o-Minimal Structures*.
London Mathematical Society Lecture Note Series. Cambridge University
Press.

# Agreement Value Significativity

Agreement measures, such as Cohen’s $\kappa$(Cohen 1960) and
$\text{IA}$(Casagrande, Fabris, and Girometti 2020), gauge the agreement
between two discrete classifiers mapping their confusion matrices into
real values.

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

Let $\mathcal{M}_{n,m}$ be the set of $n \times n$-confusion matrices
built over a dataset having size $m$, i.e., the elements of any matrix
$\mathcal{M}_{n,m}$ sum up to $m$. For any agreement measure $\sigma$,
the $\sigma$-significativity of $c$ in $\mathcal{M}_{n,m}$ is \$\$
\varrho\_{\sigma, n, m}(c) \stackrel{\tiny\textrm{def}}{=}
\frac{\left\|\left\\M \in \mathcal{M}\_{n,m} \\ \|\\ \sigma(M) \<
c\right\\\right\|}{\left\|\mathcal{M}\_{n,m}\right\|}. \$\$ The
$\sigma$-significativity is the \[0,1\]-normalised number of matrices in
$\mathcal{M}_{n,m}$ having a $\sigma$-value greater than $c$. From a
different point of view, it corresponds to the probability of choosing
by chance a matrix whose $\sigma$-value greater than $c$.

Since
$\left| \mathcal{M}_{n,m} \right| = \left( \frac{n^{2} + m - 1}{m} \right)$,
computing $\varrho_{\sigma,n,m}(c)$ takes time
$O\left( n^{2}\left( \frac{n^{2} + m - 1}{m} \right) \right)$.
Nevertheless, the Monte Carlo method (Metropolis and Ulam 1949) can
estimate $\varrho_{\sigma,n,m}(c)$ with $N$ samples in time
$\Theta\left( N\left( n^{2} + m \right) \right)$ with an error
proportional to $1/\sqrt{N}$ (e.g., see (Mackay 1998)).

The function
[`significativity()`](https://albertocasagrande.github.io/rSignificativity/reference/significativity.md)
lets us to both exactly compute $\varrho_{\sigma,n,m}(c)$ and estimate
it by using the Monte Carlo method.

``` r
library("rSignificativity")

# evaluate kappa-significativity of 0.5 in M_{2,100} with 10000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 100)
#> [1] 0.888

# evaluate kappa-significativity of 0.5 in M_{2,100} with 1000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 100, number_of_samples = 1000)
#> [1] 0.897

# evaluate kappa-significativity of 0.5 in M_{2,5} with 1000 samples
significativity(cohen_kappa, c = 0.5, n = 2, m = 5, number_of_samples = 1000)
#> [1] 0.835

# exactly compute kappa-significativity of 0.5 in M_{2,5}
significativity(cohen_kappa, c = 0.5, n = 2, m = 5, number_of_samples = NULL)
#> [1] 0.8214286
```

The $\sigma$-significativity of $c$ in $\mathcal{P}_{n}$, where
$\mathcal{P}_{n}$ is the set of all the $n \times n$-probability
matrices, is \$\$ \rho\_{\sigma, n}(c) \stackrel{\tiny\textrm{def}}{=}
\frac{V\left(\left\\M \in \mathcal{P}\_{n} \\ \|\\ \sigma(M) \<
c\right\\\right)}{V(\Delta^{(n^2-1)})} \$\$ where
$\Delta^{(k - 1)} = \left\{ \langle x_{1},\ldots x_{k}\rangle \in {\mathbb{R}}_{\geq 0}^{k}\,|\,\sum_{i = 1}^{k}x_{i} = 1 \right\}$
is the $k$-dimensional probability simplex and $V( \cdot )$ is the
$n^{2} - 1$-dimensional Lebesgue measure (Casagrande et al. 2025).

If $\sigma(M) < c$ is definable in an o-minimal theory (van den Dries
1998), such as in the case of Cohen’s $\kappa$ and $\text{IA}*$, then
$$\lim\limits_{m\rightarrow + \infty}\varrho_{\sigma,n,m}(c) = \rho_{\sigma,n}(c).$$

We can estimate $\rho_{\sigma,n}(c)$ by using the Monte Carlo method
with $N$ samples in time $\Theta\left( n^{2}N \right)$. The function
[`significativity()`](https://albertocasagrande.github.io/rSignificativity/reference/significativity.md)
also implements this algorithm.

``` r
# evaluate kappa-significativity of 0.5 in P_{2} with 1000 samples
significativity(cohen_kappa, 0.5, n = 2, number_of_samples = 1000)
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

Casagrande, Alberto, Francesco Fabris, and Rossano Girometti. 2020.
“Extending Information Agreement by Continuity.” In *2020 IEEE
International Conference on Bioinformatics and Biomedicine (BIBM)*,
1432–39. IEEE. <https://doi.org/10.1109/BIBM49941.2020.9313173>.

Casagrande, Alberto, Francesco Fabris, Girometti Rossano, and Roberto
Pagliarini. 2025. “Statical Coefficient Significativity.”

Cohen, Jacob. 1960. “A Coefficient of Agreement for Nominal Scales.”
*Educational and Psychological Measurement* 20: 37–46.
<https://doi.org/10.1177/001316446002000104>.

Mackay, D. J. C. 1998. “Introduction to Monte Carlo Methods.” In
*Learning in Graphical Models*, edited by Michael I. Jordan, 175–204.
Dordrecht: Springer Netherlands.
<https://doi.org/10.1007/978-94-011-5014-9_7>.

Metropolis, Nicholas, and Stanisław Ulam. 1949. “The Monte Carlo
Method.” *Journal of the American Statistical Association* 44 (247):
335–41. <https://doi.org/10.1080/01621459.1949.10483310>.

van den Dries, Lou P. D. 1998. *Tame Topology and o-Minimal Structures*.
London Mathematical Society Lecture Note Series. Cambridge University
Press.

# Cohen's \\\kappa\\

This method computes Cohen's \\\kappa\\ of a confusion matrix.

## Usage

``` r
cohen_kappa(conf_matrix)
```

## Arguments

- conf_matrix:

  The confusion matrix.

## Value

The Cohen's \\\kappa\\ of `conf_matrix`.

## Examples

``` r
# create a confusion matrix
conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)

# evaluate Cohen's kappa of `conf_matrix`
cohen_kappa(conf_matrix)
#> [1] -0.04166667
```

# Scott's \\\pi\\

This method computes Scott's \\\pi\\ of a confusion matrix.

## Usage

``` r
scott_pi(conf_matrix)
```

## Arguments

- conf_matrix:

  The confusion matrix.

## Value

The Scott's \\\pi\\ of `conf_matrix`.

## Examples

``` r
# create a confusion matrix
conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)

# evaluate Scott's pi of `conf_matrix`
scott_pi(conf_matrix)
#> [1] -0.05633803
```

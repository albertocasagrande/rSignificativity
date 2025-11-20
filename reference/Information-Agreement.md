# Information Agreement

This method computes the Information Agreement \\\text{IA}\\ of a
confusion matrix.

## Usage

``` r
IA(conf_matrix)
```

## Arguments

- conf_matrix:

  The confusion matrix.

## Value

The \\\text{IA}\\ of `conf_matrix`.

## Examples

``` r
# create a confusion matrix
conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)

# evaluate the Information Agreement of `conf_matrix`
IA(conf_matrix)
#> [1] 0.005631984
```

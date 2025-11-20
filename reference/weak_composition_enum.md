# The lexicographic weak composition enumeration

This function implements the lexicographic order enumerator for weak
compositions of \\m\\ in \\k\\ parts.

A weak composition of \\m\\ in \\k\\ parts is a vector of \\k\\ natural
numbers whose sum is \\m\\. \\\mathcal{C}\_{m,k}\\ is the set of all the
weak composition of \\m\\ in \\k\\ parts.

The lexicographic order \\\<\_l\\ over \\\mathcal{C}\_{m,k}\\ is an
order between weak compositions of \\m\\ in \\k\\ parts such that
\\\langle a_1,\ldots a_k \rangle \<\_l \langle b_1,\ldots b_k \rangle\\
when there exist \\i \in \[1, k\]\\ such that for all \\j \in \[1,
i-1\]\\ it holds that \\a_j = b_j\\ and \\a_i\<b_i\\.

## Arguments

- m:

  The sum of the natural values in the output vector.

- k:

  The size of the output vector.

- index:

  The index of the output vector in the lexicographic order.

## Value

The `index`-th weak composition in the lexicographic order of
\\\mathcal{C}\_{m,k}\\.

## Examples

``` r
# evaluate the 1st element in the lexicographic order of C_{100, 5}
weak_composition_enum(100, 5, 0)
#> [1]   0   0   0   0 100

# evaluate it again
weak_composition_enum(100, 5, 0)
#> [1]   0   0   0   0 100

# evaluate the 42137th element in the lexicographic order of C_{100, 5}
weak_composition_enum(100, 5, 42137)
#> [1]  0  8 56 33  3
```

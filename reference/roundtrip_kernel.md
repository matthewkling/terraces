# Round-trip kernel for a disaggregation method

For prefilter methods, returns the analytical kernel representing
"disagg then aggregate" as a small convolution on the coarse grid.
Useful for diagnostics and for understanding the bias structure of
standard interpolation.

## Usage

``` r
roundtrip_kernel(method, fact)
```

## Arguments

- method:

  Character. Currently \`"bilinear"\` (3x3 kernel) and \`"cubic"\` (5x5
  kernel, Keys convention with a = -0.5) are supported as prefilter
  methods. (Pycnophylactic disaggregation is iterative and does not have
  a round-trip kernel in this sense.)

- fact:

  Integer disagg factor.

## Value

A square odd-sized numeric matrix.

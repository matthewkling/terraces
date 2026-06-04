# Compute the inverse kernel for a prefilter disaggregation method

The inverse kernel is the focal filter applied to the coarse raster
before standard interpolation; it predistorts the coarse data so that
subsequent disagg produces a fine raster whose block means equal the
input coarse values.

## Usage

``` r
kernel(method, fact, radius = NULL)
```

## Arguments

- method:

  Character, prefilter method name. Currently \`"bilinear"\`.

- fact:

  Integer disagg factor.

- radius:

  Integer half-width of the kernel. Larger = more accurate at higher
  one-time cost. \`NULL\` uses a method-specific default.

## Value

A square matrix of side \`2 \* radius + 1\`, with attributes \`method\`,
\`fact\`, \`radius\`, \`tail_max\` (largest absolute value at the kernel
boundary; an upper bound on truncation error), and \`conv_max_err\`
(largest deviation of K \* K_inv from the delta on the interior; should
be machine epsilon).

## Details

Kernels depend only on the method and disagg factor (not the data) and
are cached for the R session.

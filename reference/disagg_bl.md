# Pre-sharpened bilinear disaggregation

Disaggregates \`coarse\` to a finer resolution by integer factor
\`fact\`, such that aggregating the result back by block mean recovers
the input coarse values up to a small truncation error governed by
\`radius\`. The fine raster is a true bilinear surface — same smoothness
as standard bilinear disaggregation, but mass-preserving rather than
treating coarse values as point samples at cell centers.

## Usage

``` r
disagg_bl(
  coarse,
  fact,
  radius = NULL,
  max_radius_frac = 1/3,
  na_fill = c("auto", "reflect", "fill")
)
```

## Arguments

- coarse:

  SpatRaster. Multi-layer rasters are supported (the kernel is applied
  to each layer in turn).

- fact:

  Integer disagg factor (\>= 2).

- radius:

  Integer half-width of the inverse kernel. If \`NULL\` (the default),
  \`radius\` is set to \`max(7, fact + 2)\`, which targets a relative
  round-trip error of roughly 1e-5 to 1e-7 in the interior. Users rarely
  need to change this. Larger values reduce truncation error further at
  a small one-time cost (the kernel is cached per \`fact\`); smaller
  values are faster on small rasters but allow larger residual error.
  Automatically reduced for small rasters; see \`max_radius_frac\`.

- max_radius_frac:

  Numeric. Upper bound on radius as a fraction of the coarse raster's
  smaller dimension. Defaults to \`1/3\` to keep the boundary band
  (where reflective extension is approximate) from dominating the
  interior. Allowed range: (0, 0.5\].

- na_fill:

  Boundary/NA handling mode passed to \[apply_kernel()\]. \`"auto"\`
  (default) selects \`"reflect"\` for rasters with no NAs and \`"fill"\`
  for rasters with NAs.

## Value

Fine SpatRaster. Layer names preserved.

## Details

Internally: 1. Computes (or looks up cached) the bilinear inverse kernel
for \`fact\`. 2. Applies the kernel as a focal pass on \`coarse\`,
pre-sharpening the coarse raster so that the subsequent standard
bilinear disaggregation produces a result whose block means equal the
input. 3. Bilinear-disaggregates the pre-sharpened coarse raster via
\[terra::disagg()\].

\*\*Mass preservation is approximate.\*\* The inverse kernel is a finite
approximation to an ideal infinite operator, so a small truncation error
remains. The error is concentrated in cells within \`radius\` of the
raster boundary; interior cells have negligible error. Default radii
target ~1e-5 to ~1e-7 relative round-trip error in the interior. Setting
a larger \`radius\` reduces the error further at a small one-time cost
(the kernel is cached per disaggregation factor).

\*\*Boundary handling.\*\* The focal pass uses reflective extension at
raster edges by default. For rasters with NAs, the \`"fill"\` mode pads
with nearest-valid values instead. See \[edge_effects()\] to visualize
which fine cells fall in the edge-affected zone for given inputs.

## References

Unser, M. (1999). Splines: A perfect fit for signal and image
processing. \*IEEE Signal Processing Magazine\* 16(6), 22-38.

## Examples

``` r
library(terra)
#> terra 1.9.27
coarse <- rast(nrows = 30, ncols = 30, vals = runif(900))

# standard bilinear: not mass-preserving
fine_std <- disagg(coarse, fact = 5, method = "bilinear")
back_std <- aggregate(fine_std, fact = 5, fun = "mean")
max(abs(values(coarse) - values(back_std)))   # substantial
#> [1] 0.3261633

# pre-sharpened bilinear: mass-preserving up to truncation error
fine_bl  <- disagg_bl(coarse, fact = 5)
back_bl  <- aggregate(fine_bl, fact = 5, fun = "mean")
max(abs(values(coarse) - values(back_bl)))    # small; reduced by larger radius
#> [1] 0.1625402
```

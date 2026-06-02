#' Pre-sharpened bilinear disaggregation
#'
#' Disaggregates `coarse` to a finer resolution by integer factor `fact`,
#' such that aggregating the result back by block mean recovers the input
#' coarse values up to a small truncation error governed by `radius`. The
#' fine raster is a true bilinear surface — same smoothness as standard
#' bilinear disaggregation, but mass-preserving rather than treating coarse
#' values as point samples at cell centers.
#'
#' Internally:
#' 1. Computes (or looks up cached) the bilinear inverse kernel for `fact`.
#' 2. Applies the kernel as a focal pass on `coarse`, pre-sharpening the
#'    coarse raster so that the subsequent standard bilinear disaggregation
#'    produces a result whose block means equal the input.
#' 3. Bilinear-disaggregates the pre-sharpened coarse raster via
#'    [terra::disagg()].
#'
#' @param coarse SpatRaster. Multi-layer rasters are supported (the kernel
#'   is applied to each layer in turn).
#' @param fact Integer disagg factor (>= 2).
#' @param radius Integer half-width of the inverse kernel. If `NULL`
#'   (the default), `radius` is set to `max(7, fact + 2)`, which targets
#'   a relative round-trip error of roughly 1e-5 to 1e-7 in the interior.
#'   Users rarely need to change this. Larger values reduce truncation
#'   error further at a small one-time cost (the kernel is cached per
#'   `fact`); smaller values are faster on small rasters but allow larger
#'   residual error. Automatically reduced for small rasters; see
#'   `max_radius_frac`.
#' @param max_radius_frac Numeric. Upper bound on radius as a fraction of
#'   the coarse raster's smaller dimension. Defaults to `1/3` to keep the
#'   boundary band (where reflective extension is approximate) from
#'   dominating the interior. Allowed range: (0, 0.5].
#' @param na_fill Boundary/NA handling mode passed to [apply_kernel()].
#'   `"auto"` (default) selects `"reflect"` for rasters with no NAs and
#'   `"fill"` for rasters with NAs.
#'
#' @return Fine SpatRaster. Layer names preserved.
#'
#' @details
#' **Mass preservation is approximate.** The inverse kernel is a finite
#' approximation to an ideal infinite operator, so a small truncation
#' error remains. The error is concentrated in cells within `radius` of
#' the raster boundary; interior cells have negligible error. Default
#' radii target ~1e-5 to ~1e-7 relative round-trip error in the interior.
#' Setting a larger `radius` reduces the error further at a small one-time
#' cost (the kernel is cached per disaggregation factor).
#'
#' **Boundary handling.** The focal pass uses reflective extension at
#' raster edges by default. For rasters with NAs, the `"fill"` mode pads
#' with nearest-valid values instead. See [edge_effects()] to visualize
#' which fine cells fall in the edge-affected zone for given inputs.
#'
#' @references
#' Unser, M. (1999). Splines: A perfect fit for signal and image processing.
#' *IEEE Signal Processing Magazine* 16(6), 22-38.
#'
#' @examples
#' library(terra)
#' coarse <- rast(nrows = 30, ncols = 30, vals = runif(900))
#'
#' # standard bilinear: not mass-preserving
#' fine_std <- disagg(coarse, fact = 5, method = "bilinear")
#' back_std <- aggregate(fine_std, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back_std)))   # substantial
#'
#' # pre-sharpened bilinear: mass-preserving up to truncation error
#' fine_bl  <- disagg_bl(coarse, fact = 5)
#' back_bl  <- aggregate(fine_bl, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back_bl)))    # small; reduced by larger radius
#'
#' @export
disagg_bl <- function(coarse, fact, radius = NULL,
                      max_radius_frac = 1/3,
                      na_fill = c("auto", "reflect", "fill")) {
      na_fill <- match.arg(na_fill)
      .validate_inputs(coarse, fact, max_radius_frac)
      fact <- as.integer(fact)

      if (is.null(radius)) {
            radius <- .default_radius("bilinear", fact)
      }
      radius <- as.integer(radius)

      min_dim <- min(nrow(coarse), ncol(coarse))
      r_hard <- (min_dim - 1L) %/% 2L
      r_soft <- floor(max_radius_frac * min_dim)
      r_max  <- min(r_hard, r_soft)
      if (r_max < 2L) {
            stop("Coarse raster too small (", nrow(coarse), " x ", ncol(coarse),
                 ") for this method; need at least ~6 cells per side.")
      }
      if (radius > r_max) {
            warning("Inverse-kernel radius reduced from ", radius, " to ", r_max,
                    " to fit coarse raster (", nrow(coarse), " x ", ncol(coarse),
                    "); round-trip accuracy will be lower than requested.")
            radius <- r_max
      }

      K_inv <- kernel("bilinear", fact, radius)
      apply_kernel(coarse, K_inv, fact, method = "bilinear", na_fill = na_fill)
}
